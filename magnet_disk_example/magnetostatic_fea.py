from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from petsc4py import PETSc
from msh2xdmf import import_mesh
from magnet_disk_mesh import generateMeshMovement, getInitialEdgeCoords
from piecewise_permeability import *
from scipy.spatial import KDTree

#set_log_level(1)
def m2p(A):
    return as_backend_type(A).mat()

def v2p(v):
    return as_backend_type(v).vec()

def transpose(A):
    """
    Transpose for matrix of DOLFIN type
    """
    return PETScMatrix(as_backend_type(A).mat().transpose(PETSc.Mat(MPI.comm_world)))

def computeMatVecProductFwd(A, x):
    """
    Compute y = A * x
    A: ufl form matrix
    x: ufl function
    """
    A_p = m2p(A)
    y = A_p * v2p(x.vector())
    y.assemble()
#    y.ghostUpdate()
    return y.getArray()


def computeMatVecProductBwd(A, R):
    """
    Compute y = A.T * R
    A: ufl form matrix
    R: ufl function
    """
    row, col = m2p(A).getSizes()
    y = PETSc.Vec().create()
    y.setSizes(col)
    y.setUp()
    m2p(A).multTranspose(v2p(R.vector()),y)
    y.assemble()
#    y.ghostUpdate()
    return y.getArray()

def update(f, f_values):
    """
    f: dolfin function
    f_values: numpy array
    """
    f.vector().set_local(f_values)
    v2p(f.vector()).assemble()
    v2p(f.vector()).ghostUpdate()
    
def findNodeIndices(node_coordinates, coordinates):
    tree = KDTree(coordinates)
    dist, node_indices = tree.query(node_coordinates)
    return node_indices
    
I = Identity(2)
def gradx(f,u):
    return dot(grad(f), inv(I + grad(u)))

def J(uhat):
    return det(I + grad(uhat))
    
# START NEW PERMEABILITY
def RelativePermeability(subdomain, u, uhat):
    gradu = gradx(u,uhat)
    if subdomain == 1: # Electrical/Silicon/Laminated Steel
        B = as_vector((gradu[1], -gradu[0]))
        norm_B = sqrt(dot(B, B) + DOLFIN_EPS)
        
        mu = conditional(
            lt(norm_B, 1.004),
            linearPortion(norm_B),
            conditional(
                lt(norm_B, 1.433),
                cubicPortion(norm_B),
                (expA * exp(expB*norm_B + expC) + 1)
            )
        )
        # mu = 4000. # constant case

    else: # Magnets
        mu = 1.05

    return mu

# END NEW PERMEABILITY

def Jm(v,uhat,num_magnets,Hc):
    Jm = 0.
    gradv = gradx(v,uhat)
    for i in range(num_magnets):
        angle = np.pi / 4 + i * np.pi / 2
        Hx = Constant((-1)**i * Hc * np.cos(angle))
        Hy = Constant((-1)**i * Hc * np.sin(angle))

        H = as_vector([Hx, Hy])
        
        curl_v = as_vector([gradv[1],-gradv[0]])
        Jm += inner(H,curl_v)*dx(i+2)
    return Jm
    
def pdeRes(u,v,uhat,dx,num_magnets,Hc,vacuum_perm):
    res = 0
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    for i in range(5):
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *dot(gradu,gradv)*dx(i + 1)
    res -= Jm(v,uhat,num_magnets,Hc)
    return res


class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
        

class MagnetostaticProblem(object):
    """
    Preprocessor to set up the mesh and mesh functions
    """  
    def __init__(self, mesh_file="magnet_disk_mesh"):

        self.mesh_file = mesh_file
        self.initMesh()
        self.initFunctionSpace()
        self.locateMeshMotion()
        
        # PROBLEM SPECIFIC PARAMETERS
        self.Hc = 838.e3  # 838 kA/m
        self.num_magnets = 4
        self.vacuum_perm = 4e-7 * np.pi
        
        self.uhat = Function(self.VHAT) # Function for the solution of the mesh motion
        self.uhat_old = Function(self.VHAT) # Function to apply the BCs of the mesh motion
        self.A_z = Function(self.V) # Function for the solution of the magnetic vector potential
        self.v = TestFunction(self.V)
        self.dR = Function(self.V)
        self.du = Function(self.V)
        self.duhat = Function(self.VHAT)
        
        self.dR_du = derivative(self.resMS(), self.A_z)
        self.dR_duhat = derivative(self.resMS(), self.uhat)
        
        self.total_dofs_A_z = len(self.A_z.vector().get_local())
        self.total_dofs_uhat = len(self.uhat.vector().get_local())
        
    def initMesh(self):
        self.mesh, self.boundaries_mf, self.subdomains_mf, association_table = import_mesh(
            prefix=self.mesh_file,
            dim=2,
            subdomains=True
        )
        self.dx = Measure('dx', domain=self.mesh, subdomain_data=self.subdomains_mf)
        self.dS = Measure('dS', domain=self.mesh, subdomain_data=self.boundaries_mf)
        
    def initFunctionSpace(self):
        self.V = FunctionSpace(self.mesh, 'P', 1)
        self.VHAT = VectorFunctionSpace(self.mesh, 'P', 1)

    
    def locateMeshMotion(self):
        V0 = FunctionSpace(self.mesh, 'CG', 1)
        coordinates = V0.tabulate_dof_coordinates()
        
        old_edge_coords = getInitialEdgeCoords()
        
        node_indices = findNodeIndices(np.reshape(old_edge_coords, (-1,2)), 
                                        coordinates)
        edge_indices = np.empty(2*len(node_indices))
        for i in range(len(node_indices)):
            edge_indices[2*i] = 2*node_indices[i]
            edge_indices[2*i+1] = 2*node_indices[i]+1
        self.edge_indices = edge_indices
    
    def resM(self):
        # Residual for mesh, which satisfies a fictitious elastic problem:
        duhat = TestFunction(self.VHAT)
        F_m = grad(self.uhat) + I
        E_m = 0.5*(F_m.T*F_m - I)
        m_jac_stiff_pow = Constant(3)
        # Artificially stiffen the mesh where it is getting crushed:
        K_m = Constant(1)/pow(det(F_m),m_jac_stiff_pow)
        mu_m = Constant(1)/pow(det(F_m),m_jac_stiff_pow)
        S_m = K_m*tr(E_m)*I + 2.0*mu_m*(E_m - tr(E_m)*I/3.0)
        res_m = (inner(F_m*S_m,grad(duhat)))*dx
        return res_m
        
    def resMS(self):
        res_ms = pdeRes(
                self.A_z,self.v,self.uhat,self.dx,
                self.num_magnets,self.Hc,self.vacuum_perm)
        return res_ms
        
    def setBCMagnetostatic(self):
        A_outer = Constant(0.0)
        outer_bound = OuterBoundary()
        return DirichletBC(self.V, self.A_z, outer_bound)
    
    def setBCMeshMotion(self, edge_deltas):
        return DirichletBC(self.VHAT, self.uhat_old, self.boundaries_mf, 1000)

    def getDisplacementSteps(self, edge_deltas):
        STEPS = 2
        max_disp = np.max(np.abs(edge_deltas))
        min_cell_size = self.mesh.hmin()
        min_STEPS = round(max_disp/min_cell_size)
        if min_STEPS >= STEPS:
            STEPS = min_STEPS
        return STEPS
        
    def solveMeshMotion(self, edge_deltas=None, STEPS=2):
        if edge_deltas is None:
            edge_deltas = generateMeshMovement()

        STEPS = self.getDisplacementSteps(edge_deltas)
        self.increment_deltas = edge_deltas/STEPS
        self.edge_deltas = edge_deltas
        
        ####### Formulation of mesh motion as a hyperelastic problem #######
        
        # Initialize the boundary condition for the first displacement step
        self.uhat_old.vector()[self.edge_indices] = self.increment_deltas
        bc_m = self.setBCMeshMotion(self.uhat_old)
        res_m = self.resM()
        Dres_m = derivative(res_m, self.uhat)

        ####### Nonlinear solver setup #######

        # Nonlinear solver parameters
        REL_TOL_M = 1e-4
        MAX_ITERS_M = 100

        # Set up nonlinear problem for mesh motion:
        problem_m = NonlinearVariationalProblem(res_m, self.uhat, 
                                                bc_m, Dres_m)
        solver_m = NonlinearVariationalSolver(problem_m)
        solver_m.parameters['newton_solver']\
                           ['maximum_iterations'] = MAX_ITERS_M
        solver_m.parameters['newton_solver']\
                           ['relative_tolerance'] = REL_TOL_M
        
        for i in range(STEPS):
            print(80*"=")
            print("  FEA: Step "+str(i+1)+" of mesh movement")
            print(80*"=")
            solver_m.solve()
            self.uhat_old.assign(self.uhat)
            self.uhat_old.vector()[self.edge_indices] += self.increment_deltas
            
            
        print(80*"=")
        print(' FEA: L2 error of the mesh motion on the edges:', 
                    np.linalg.norm(self.uhat.vector()[self.edge_indices]
                                     - self.edge_deltas))
        print(80*"=")
        
    def solveMagnetostatic(self):
        
        ######### Formulation of the magnetostatic problem #################
        
        print(80*"=")
        print(" FEA: Solving the magnetostatic problem on the deformed mesh")
        print(80*"=")
        bc_ms = self.setBCMagnetostatic()
        solve(self.resMS()==0, self.A_z, J=self.dR_du, bcs=bc_ms)
        self.B = project(as_vector((self.A_z.dx(1),-self.A_z.dx(0))),
                        VectorFunctionSpace(self.mesh,'DG',0))

    def solveLinearFwd(self, A, dR):
        """
        solve linear system dR = dR_du (A) * du
        """
        self.dR.vector().set_local(dR)
        v2p(self.dR.vector()).assemble()
        v2p(self.dR.vector()).ghostUpdate()

        self.du.vector().set_local(np.zeros(self.local_dof_u))
        v2p(self.du.vector()).assemble()
        v2p(self.du.vector()).ghostUpdate()
    
        solverFwd = LUSolver("mumps")
        solverFwd.solve(A, self.du.vector(), self.dR.vector())
        v2p(self.du.vector()).assemble()
        v2p(self.du.vector()).ghostUpdate()
        return self.du.vector().get_local()

    def solveLinearBwd(self, A, du):
        """
        solve linear system du = dR_du.T (A_T) * dR
        """
        self.du.vector().set_local(du)
        v2p(self.du.vector()).assemble()
        v2p(self.du.vector()).ghostUpdate()

        self.dR.vector().set_local(np.zeros(self.local_dof_u))
        v2p(self.dR.vector()).assemble()
        v2p(self.dR.vector()).ghostUpdate()
        
        A_T = transpose(A)

        solverBwd = LUSolver("mumps")
        solverBwd.solve(A_T, self.dR.vector(), self.du.vector())
        v2p(self.dR.vector()).assemble()
        v2p(self.dR.vector()).ghostUpdate()
        return self.dR.vector().get_local()

    def moveMesh(self):
        ALE.move(self.mesh, self.uhat)

if __name__ == "__main__":
    problem = MagnetostaticProblem()
    problem.solveMeshMotion()
#    problem.solveMagnetostatic(edge_deltas=np.zeros(len(problem.edge_indices)))
#    problem.solveMagnetostatic()
    plt.figure(1)
    #plot(problem.B, linewidth=40)
    ALE.move(problem.mesh, problem.uhat)
    plot(problem.mesh)
#    print(problem.A_z.vector().get_local()[:10])
    plt.show()




