from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from petsc4py import PETSc
from msh2xdmf import import_mesh
# from magnet_disk_mesh import generateMeshMovement, getInitialEdgeCoords
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
    # if subdomain == 1: # Electrical/Silicon/Laminated Steel
    if subdomain == 1 or subdomain == 2: # Electrical/Silicon/Laminated Steel
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
    elif subdomain == 3: 
        mu = 1.00 # insert value for titanium or shaft material
    elif subdomain >= 4 and subdomain <= 28: # AIR
        mu = 1.0
    elif subdomain >= 29 and subdomain <= 40: # NEODYMIUM
        mu = 1.05
    elif subdomain >= 41: # COPPER
        mu = 1.00

    return mu

# END NEW PERMEABILITY

def JS(v,uhat,p,s,Hc):
    Jm = 0.
    gradv = gradx(v,uhat)
    base_magnet_dir = 2 * np.pi / p / 2
    magnet_sweep    = 2 * np.pi / p
    for i in range(p):
        angle = base_magnet_dir + i * magnet_sweep
        Hx = Constant((-1)**i * Hc * np.cos(angle))
        Hy = Constant((-1)**i * Hc * np.sin(angle))

        H = as_vector([Hx, Hy])
        
        curl_v = as_vector([gradv[1],-gradv[0]])
        Jm += inner(H,curl_v)*dx(i + 4 + p*2 + 1)

    num_windings = 2*s
    num_phases = 3
    coil_per_phase = 6
    stator_winding_index_start  = 4 + 3 * p + 1
    stator_winding_index_end    = stator_winding_index_start + num_windings
    Jw = 0.
    # N = 39 / 1e-3
    N = 39
    JA, JB, JC = 94.26 * N, 47.13 * N, -47.13 * N
    # JA, JB, JC = 66.65 * N, -91.04 * N, 24.4 * N
    for i in range(int((num_windings) / (num_phases * coil_per_phase))):
        coil_start_ind = i * num_phases * coil_per_phase
        for j in range(3):
            if j == 0:
                J = JA
            elif j == 1:
                J = JB
            elif j == 2:
                J = JC
            phase_start_ind = coil_per_phase * j
            J_list = [
                J * (-1)**i * v * dx(stator_winding_index_start + phase_start_ind + coil_start_ind),
                J * (-1)**(i+1) * v * dx(stator_winding_index_start + 1 + phase_start_ind + coil_start_ind),
                J * (-1)**(i+1) * v * dx(stator_winding_index_start + 2 + phase_start_ind + coil_start_ind),
                J * (-1)**i * v * dx(stator_winding_index_start + 3 + phase_start_ind + coil_start_ind),
                J * (-1)**i * v * dx(stator_winding_index_start + 4 + phase_start_ind + coil_start_ind),
                J * (-1)**(i+1) * v * dx(stator_winding_index_start + 5 + phase_start_ind + coil_start_ind),
            ]
            Jw += sum(J_list)
            # order: + - - + + - (signs switch with each instance of the phases)
    
    
    return Jm + Jw
        
def pdeRes(u,v,uhat,dx,p,s,Hc,vacuum_perm):
    res = 0.
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    num_components = 4 * 3 * p + 2 * s
    for i in range(num_components):
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *dot(gradu,gradv)*dx(i + 1)
    res -= JS(v,uhat,p,s,Hc)
    return res


class OuterBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
        
class MagnetostaticProblem(object):
    """
    Preprocessor to set up the mesh and mesh functions
    """  
    def __init__(self, mesh_file="motor_mesh_1"):

        self.mesh_file = mesh_file
        self.initMesh()
        self.initFunctionSpace()
        # self.locateMeshMotion()
        
        # PROBLEM SPECIFIC PARAMETERS
        self.Hc = 838.e3  # 838 kA/m
        self.p = 12
        self.s = 3 * self.p
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
                self.p,self.s,self.Hc,self.vacuum_perm)
        return res_ms

    def getSubdomainArea(self, func_unit, subdomain):
        area = assemble(inner(func_unit, func_unit)*self.dx((subdomain)))
        return area

    def getFuncAverageSubdomain(self, func, subdomain):
        func_unit = interpolate(Constant(1.0), func.function_space())
        integral = assemble(inner(func, func_unit)*self.dx((subdomain)))
        area = self.getSubdomainArea(func_unit, subdomain)
        # area = assemble(inner(func_unit, func_unit)*self.dx(subdomain))
        avg_func = integral/area
        return avg_func
     
     # TODO: add the formula of flux linkage using 'getFuncAverageSubdomain' for each winding
    def extractSubdomainAverageA_z(self, func, subdomain_range):
        subdomain_avg_A_z = []
        subdomain_avg_A_z_deltas = []
        # winding_start, winding_end = 41, 112
        for subdomain in subdomain_range:
            subdomain_avg_A_z.append(
                abs(self.getFuncAverageSubdomain(func, subdomain))
            )
        for i in range(int(len(subdomain_range)/2)):
            subdomain_avg_A_z_deltas.append(
                abs(subdomain_avg_A_z[2*i] - subdomain_avg_A_z[2*i+1])
            )
        
        return subdomain_avg_A_z, subdomain_avg_A_z_deltas

    def fluxLinkage(self):
        pass
        
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
        res_ms = self.resMS()
        A_z_ = TrialFunction(self.V)
        Dres_ms = derivative(res_ms, self.A_z, A_z_)
        
       
        # Nonlinear solver parameters
        REL_TOL_M = 1e-4
        MAX_ITERS_M = 100

#        update(self.A_z, 0.1*np.ones(self.total_dofs_A_z))
        problem_ms = NonlinearVariationalProblem(res_ms, self.A_z, 
                                                bc_ms, Dres_ms)
        solver_ms = NonlinearVariationalSolver(problem_ms)
        solver_ms.parameters['nonlinear_solver']='snes' 
        solver_ms.parameters['snes_solver']['line_search'] = 'bt' 
        solver_ms.parameters['snes_solver']['absolute_tolerance'] = REL_TOL_M
        solver_ms.parameters['snes_solver']['maximum_iterations'] = MAX_ITERS_M
        solver_ms.parameters['snes_solver']['linear_solver']='mumps' # "cg" "gmres"
        solver_ms.parameters['snes_solver']['error_on_nonconvergence'] = False
        solver_ms.solve()
        self.B = project(as_vector((self.A_z.dx(1),-self.A_z.dx(0))),
                        VectorFunctionSpace(self.mesh,'DG',0))

        subdomain_range = np.arange(41,112+1)
        self.winding_delta_A_z = np.array(
            self.extractSubdomainAverageA_z(
                func=self.A_z,
                subdomain_range=subdomain_range
        )[1])

        area_func_unit = interpolate(Constant(1.0), self.A_z.function_space())
        self.winding_area = self.getSubdomainArea(
            func_unit=area_func_unit,
            subdomain=42
        )

        # print(self.winding_delta_A_z)
        

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
    # problem.solveMeshMotion()
    problem.solveMagnetostatic()

    vtkfile_A_z = File('solutions/Magnetic_Vector_Potential.pvd')
    vtkfile_B = File('solutions/Magnetic_Flux_Density.pvd')
    vtkfile_A_z << problem.A_z
    vtkfile_B << problem.B

    ###### Test the average calculation for the flux linkage
    subdomain_range = np.arange(41,112)
    asdf, deltas  = problem.extractSubdomainAverageA_z(
        func=problem.A_z,
        subdomain_range=subdomain_range
    )
    for i in range(len(subdomain_range)):
        print("Average A_z for subdomain "+str(i+41))
        # print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
        print(asdf[i])

    for i in range(len(deltas)):
        print("Delta A_z for Stator Tooth "+str(i+1))
        # print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
        print(deltas[i])
    
    plt.figure(1)
    plot(problem.A_z)
    # plot(problem.B, linewidth=40)
#    ALE.move(problem.mesh, problem.uhat)
    # plot(problem.mesh)
#    print(problem.A_z.vector().get_local()[:10])
    plt.show()




