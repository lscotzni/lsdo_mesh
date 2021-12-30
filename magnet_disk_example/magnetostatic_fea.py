"""
The FEniCS wrapper for the variational forms and the partial derivatives computation
in the magnetostatic/motor problem on a deformable mesh.
"""


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from petsc4py import PETSc
from msh2xdmf import import_mesh
from magnet_disk_mesh import generateMeshMovement, getInitialEdgeCoords
from piecewise_permeability import *
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix


#set_log_level(1)
def m2p(A):
    """
    Convert the matrix of DOLFIN type to a PETSc.Mat object
    """
    return as_backend_type(A).mat()

def v2p(v):
    """
    Convert the vector of DOLFIN type to a PETSc.Vec object
    """
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

def zero_petsc_vec(size, comm=MPI.comm_world):
    """
    Create zero PETSc vector of size ``num_el``.
    Parameters
    ----------
    size : int
    vec_type : str, optional
        For petsc4py.PETSc.Vec types, see petsc4py.PETSc.Vec.Type.
    comm : MPI communicator
    Returns
    -------
    v : petsc4py.PETSc.Vec
    """
    v = PETSc.Vec().create(comm)
    v.setSizes(size)
    v.setUp()
    v.assemble()
    return v

def zero_petsc_mat(row, col, comm=MPI.comm_world):
    """
    Create zeros PETSc matrix with shape (``row``, ``col``).
    Parameters
    ----------
    row : int
    col : int
    mat_type : str, optional
        For petsc4py.PETSc.Mat types, see petsc4py.PETSc.Mat.Type
    comm : MPI communicator
    Returns
    -------
    A : petsc4py.PETSc.Mat
    """
    A = PETSc.Mat(comm)
    A.createAIJ([row, col], comm=comm)
    A.setUp()
    A.assemble()
    return A

def convertToDense(A_petsc):
    """
    Convert the PETSc matrix to a dense numpy array
    (super unefficient, only used for debugging purposes)
    """
    A_petsc.assemble()
    A_dense = A_petsc.convert("dense")
    return A_dense.getDenseArray()

def update(f, f_values):
    """
    Update the nodal values in every dof of the DOLFIN function `f`
    according to `f_values`.
    -------------------------
    f: dolfin function
    f_values: numpy array
    """
    f.vector().set_local(f_values)
    v2p(f.vector()).assemble()
    v2p(f.vector()).ghostUpdate()

def findNodeIndices(node_coordinates, coordinates):
    """
    Find the indices of the closest nodes, given the `node_coordinates`
    for a set of nodes and the `coordinates` for all of the vertices
    in the mesh, by using scipy.spatial.KDTree
    """
    tree = KDTree(coordinates)
    dist, node_indices = tree.query(node_coordinates)
    return node_indices

I = Identity(2)
def gradx(f,uhat):
    """
    Convert the differential operation from the reference domain
    to the measure in the deformed configuration based on the mesh
    movement of `uhat`
    --------------------------
    f: DOLFIN function for the solution of the physical problem
    uhat: DOLFIN function for mesh movements
    """
    return dot(grad(f), inv(I + grad(uhat)))

def J(uhat):
    """
    Compute the determinant of the deformation gradient used in the
    integration measure of the deformed configuration wrt the the
    reference configuration.
    ---------------------------
    uhat: DOLFIN function for mesh movements
    """
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
    """
    The variational form for the source term (current) of the
    Maxwell equation
    """
    Jm = 0.
    gradv = gradx(v,uhat)
    for i in range(num_magnets):
        angle = np.pi / 4 + i * np.pi / 2
        Hx = Constant((-1)**i * Hc * np.cos(angle))
        Hy = Constant((-1)**i * Hc * np.sin(angle))

        H = as_vector([Hx, Hy])

        curl_v = as_vector([gradv[1],-gradv[0]])
        # note that the integration measure has been changed
        # to be J(uhat)*dx
        Jm += inner(H,curl_v)*J(uhat)*dx(i+2)
    return Jm

def pdeRes(u,v,uhat,dx,num_magnets,Hc,vacuum_perm):
    """
    The variational form of the PDE residual for the magnetostatic problem
    """
    res = 0
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    for i in range(5):
        # note that the integration measure has been changed to be J(uhat)*dx
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *dot(gradu,gradv)*J(uhat)*dx(i + 1)
    res -= Jm(v,uhat,num_magnets,Hc)
    return res



class OuterBoundary(SubDomain):
    """
    Define the outer boundary of the domain
    """
    def inside(self, x, on_boundary):
        return on_boundary


class MagnetostaticProblem(object):
    """
    The class of the FEniCS wrapper for the magnetostatic problem,
    with methods to compute the variational forms, partial derivatives,
    and solve the nonlinear/linear subproblems.
    """
    def __init__(self, mesh_file="magnet_disk_mesh_1"):

        self.mesh_file = mesh_file
        # Import the initial mesh from the mesh file to FEniCS
        self.initMesh()
        # Define the function spaces based on the initial mesh
        self.initFunctionSpace()
        # Get the indices of the vertices that would move during optimization
        self.edge_indices = self.locateMeshMotion()

        # PROBLEM SPECIFIC PARAMETERS
        self.Hc = 838.e3  # 838 kA/m
        self.num_magnets = 4
        self.vacuum_perm = 4e-7 * np.pi

        # Initialize the edge movements
        self.edge_deltas = None

        self.uhat = Function(self.VHAT) # Function for the solution of the mesh motion
        self.uhat_old = Function(self.VHAT) # Function to apply the BCs of the mesh motion
        self.duhat = Function(self.VHAT) # Function used in the CSDL model
        self.dRhat = Function(self.VHAT) # Function used in the CSDL model

        self.A_z = Function(self.V) # Function for the solution of the magnetic vector potential
        self.v = TestFunction(self.V)
        self.dR = Function(self.V) # Function used in the CSDL model
        self.du = Function(self.V) # Function used in the CSDL model

        self.total_dofs_bc = len(self.edge_indices)
        self.total_dofs_A_z = len(self.A_z.vector().get_local())
        self.total_dofs_uhat = len(self.uhat.vector().get_local())

        # Partial derivatives in the magnetostatic problem
        self.dR_du = derivative(self.resMS(), self.A_z)
        self.dR_duhat = derivative(self.resMS(), self.uhat)

        # Partial derivatives in the mesh motion subproblem
        self.dRm_duhat = derivative(self.resM(), self.uhat)
        self.dRm_dedge = self.getBCDerivatives()



    def initMesh(self):
        """
        Preprocessor 1 to set up the mesh and define the integration measure with the
        imported mesh functions from GMSH.
        """
        self.mesh, self.boundaries_mf, self.subdomains_mf, association_table = import_mesh(
            prefix=self.mesh_file,
            dim=2,
            subdomains=True
        )
        self.dx = Measure('dx', domain=self.mesh, subdomain_data=self.subdomains_mf)
        self.dS = Measure('dS', domain=self.mesh, subdomain_data=self.boundaries_mf)

    def initFunctionSpace(self):
        """
        Preprocessor 2 to define the function spaces for the mesh motion (VHAT)
        and the problem solution (V)
        """
        self.V = FunctionSpace(self.mesh, 'P', 1)
        self.VHAT = VectorFunctionSpace(self.mesh, 'P', 1)

    def locateMeshMotion(self):
        """
        Find the indices of the dofs for setting up the boundary condition
        in the mesh motion subproblem
        """
        V0 = FunctionSpace(self.mesh, 'CG', 1)
        coordinates = V0.tabulate_dof_coordinates()

        # One-time computation for the initial edge coordinates from
        # the code that creates the mesh file
        old_edge_coords = getInitialEdgeCoords()

        # Use KDTree to find the node indices of the points on the edge
        # in the mesh object in FEniCS
        node_indices = findNodeIndices(np.reshape(old_edge_coords, (-1,2)),
                                        coordinates)

        # Convert the node indices to edge indices, where each node has 2 dofs
        edge_indices = np.empty(2*len(node_indices))
        for i in range(len(node_indices)):
            edge_indices[2*i] = 2*node_indices[i]
            edge_indices[2*i+1] = 2*node_indices[i]+1

        return edge_indices.astype('int')

    def resM(self):
        """
        Formulation of mesh motion as a hyperelastic problem
        """
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
        """
        Formulation of the magnetostatic problem
        """
        res_ms = pdeRes(
                self.A_z,self.v,self.uhat,self.dx,
                self.num_magnets,self.Hc,self.vacuum_perm)
        return res_ms

    def getBCDerivatives(self):
        """
        Compute the derivatives of the PDE residual of the mesh motion
        subproblem wrt the BCs, which is a fixed sparse matrix with "-1"s
        on the entries corresponding to the edge indices.
        """

        row_ind = self.edge_indices
        col_ind = np.arange(self.total_dofs_bc)
        data = -1.0*np.ones(self.total_dofs_bc)
        M = csr_matrix((data, (row_ind, col_ind)),
                        shape=(self.total_dofs_uhat, self.total_dofs_bc))
        return M

    def getFuncAverageSubdomain(self, func, subdomain_ID):
        """
        Compute the average function value on a certain subdomain
        """
        func_unit = interpolate(Constant(1.0), func.function_space())
        integral = assemble(inner(func, func_unit)*self.dx(subdomain_ID))
        area = assemble(inner(func_unit, func_unit)*self.dx(subdomain_ID))
        avg_func = integral/area
        return avg_func

     # TODO: add the formula of flux linkage using 'getFuncAverageSubdomain' for each winding
    def fluxLinkage(self):
        pass

    def setBCMagnetostatic(self):
        """
        Set zero BC on the outer boundary of the domain
        """
        A_outer = Constant(0.0)
        outer_bound = OuterBoundary()
        return DirichletBC(self.V, self.A_z, outer_bound)

    def setBCMeshMotion(self):
        """
        Set BC for the mesh motion problem on the edges dynamically,
        based on 'self.uhat_old', which would be updated during the
        optimization and the displacement steps
        """
        return DirichletBC(self.VHAT, self.uhat_old, self.boundaries_mf, 1000)

    def getDisplacementSteps(self, edge_deltas):
        """
        Divide the edge movements into steps based on the current mesh size
        """
        STEPS = 2
        max_disp = np.max(np.abs(edge_deltas))
        self.moveMesh(self.uhat)
        min_cell_size = self.mesh.rmin()
        uhat_back = Function(self.VHAT)
        uhat_back.assign(-self.uhat)
        self.moveMesh(uhat_back)
        min_STEPS = round(max_disp/min_cell_size)
        print(min_STEPS)
        if min_STEPS >= STEPS:
            STEPS = min_STEPS
        increment_deltas = edge_deltas/STEPS
        return STEPS, increment_deltas

    def setUpMeshMotionSolver(self):
        """
        Set up the Newton solver for the mesh motion subproblem, compute the step size
        for the displacement increases based on the mesh size from the previous step.
        """
        if self.edge_deltas is None:
            self.edge_deltas = generateMeshMovement(pi/36)

        # Get the relative movements from the previous step
        relative_edge_deltas = self.edge_deltas - self.uhat.vector().get_local()[self.edge_indices]
        self.STEPS, self.increment_deltas = self.getDisplacementSteps(relative_edge_deltas)

        bc_m = self.setBCMeshMotion()
        res_m = self.resM()
        Dres_m = derivative(res_m, self.uhat)

        ####### Nonlinear solver setup #######

        # Nonlinear solver parameters
        REL_TOL_M = 1e-4
        MAX_ITERS_M = 100

        # Set up nonlinear problem for mesh motion:
        problem_m = NonlinearVariationalProblem(res_m, self.uhat,
                                                bc_m, Dres_m)
        self.solver_m = NonlinearVariationalSolver(problem_m)
        self.solver_m.parameters['newton_solver']\
                           ['maximum_iterations'] = MAX_ITERS_M
        self.solver_m.parameters['newton_solver']\
                           ['relative_tolerance'] = REL_TOL_M

    def solveMeshMotion(self):

        """
        Solve for the mesh motion subproblem
        """
        # Set up the Newton solver based on the boundary conditions
        self.setUpMeshMotionSolver()
        # Incrementally set the BCs to increase to `edge_deltas`
        for i in range(self.STEPS):
            print(80*"=")
            print("  FEA: Step "+str(i+1)+" of mesh movement")
            print(80*"=")
#            print(self.uhat.vector().get_local()[problem.edge_indices[:5]])
            self.uhat_old.assign(self.uhat)
            for ind in range(self.total_dofs_bc):
                self.uhat_old.vector()[self.edge_indices[ind]] += self.increment_deltas[ind]
            self.solver_m.solve()

        print(80*"=")
        print(' FEA: L2 error of the mesh motion on the edges:',
                    np.linalg.norm(self.uhat.vector()[self.edge_indices]
                                     - self.edge_deltas))
        print(80*"=")

    def solveMagnetostatic(self):

        """
        Solve the magnetostatic problem with the mesh movements `uhat`
        """

        print(80*"=")
        print(" FEA: Solving the magnetostatic problem on the deformed mesh")
        print(80*"=")
        bc_ms = self.setBCMagnetostatic()
        solve(self.resMS()==0, self.A_z, J=self.dR_du, bcs=bc_ms)
        self.B = project(as_vector((self.A_z.dx(1),-self.A_z.dx(0))),
                        VectorFunctionSpace(self.mesh,'DG',0))

    def solveLinearFwd(self, A, dR):
        """
        solve linear system dR = dR_du (A) * du in DOLFIN type
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
        solve linear system du = dR_du.T (A_T) * dR in DOLFIN type
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

    def moveMesh(self, disp):
        """
        Move the mesh object based on the `disp` function,
        for post-processing purposes.
        """
        ALE.move(self.mesh, disp)

if __name__ == "__main__":
    problem = MagnetostaticProblem()
    problem.edge_deltas = generateMeshMovement(-pi/48)
    problem.solveMeshMotion()
    problem.edge_deltas = generateMeshMovement(-pi/24)
    problem.solveMeshMotion()
    print(problem.uhat.vector().get_local()[problem.edge_indices[:5]])
    print(problem.edge_deltas[:5])
    # problem.solveMagnetostatic()
    plt.figure(1)
    problem.moveMesh(problem.uhat)
    plot(problem.mesh)
    # plot(problem.subdomains_mf)
    # plt.figure(2)
    # plot(problem.mesh,linewidth=0.4)
    # plot(problem.A_z)
    plt.show()

    # vtkfile_A_z = File('magnet_disk_solutions/Magnetic_Vector_Potential.pvd')
    # vtkfile_B = File('magnet_disk_solutions/Magnetic_Flux_Density.pvd')
    # vtkfile_A_z << problem.A_z
    # vtkfile_B << problem.B
#### Test the partial derivatives of the mesh motion subproblem
    # dRm_dedge = problem.getBCDerivatives()
    # print(dRm_dedge.dot(np.ones(problem.total_dofs_bc)))
###### Test the average calculation for the flux linkage
#    for i in range(4):
#        print("Average A_z for magnet "+str(i+1))
#        print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+2))
