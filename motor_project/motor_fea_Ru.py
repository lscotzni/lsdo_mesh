"""
The FEniCS wrapper for the variational forms and the partial derivatives computation
in the magnetostatic problem for the motor on a deformable mesh.
"""

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from petsc4py import PETSc
from msh2xdmf import import_mesh
# from magnet_disk_mesh import generateMeshMovement, getInitialEdgeCoords
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

def JS(v,uhat,p,s,Hc,i_abc):
    """
    The variational form for the source term (current) of the
    Maxwell equation
    """
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
    N = 39
    # JA, JB, JC = 94.26 * N, 47.13 * N, -47.13 * N
    # JA, JB, JC = 66.65 * N, -91.04 * N, 24.4 * N
    JA, JB, JC = i_abc[0] * N + DOLFIN_EPS, i_abc[1] * N + DOLFIN_EPS, i_abc[2] * N + DOLFIN_EPS
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

def pdeRes(u,v,uhat,dx,p,s,Hc,vacuum_perm,i_abc):
    """
    The variational form of the PDE residual for the magnetostatic problem
    """
    res = 0.
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    num_components = 4 * 3 * p + 2 * s
    for i in range(num_components):
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *dot(gradu,gradv)*dx(i + 1)
    res -= JS(v,uhat,p,s,Hc,i_abc)
    return res


class OuterBoundary(SubDomain):
    """
    Define the subdomain for the outer boundary
    """
    def inside(self, x, on_boundary):
        return on_boundary

class MotorProblem(object):
    """
    The class of the FEniCS wrapper for the motor problem,
    with methods to compute the variational forms, partial derivatives,
    and solve the nonlinear/linear subproblems.
    """
    def __init__(self, mesh_file="motor_mesh_1", i_abc=[0., 0., 0.]):

        self.mesh_file = mesh_file
        # Import the initial mesh from the mesh file to FEniCS
        self.initMesh()
        # Define the function spaces based on the initial mesh
        self.initFunctionSpace()
        # Get the indices of the vertices that would move during optimization
        # self.edge_indices = self.locateMeshMotion()

        # PROBLEM SPECIFIC PARAMETERS
        self.Hc = 838.e3  # 838 kA/m
        self.p = 12
        self.s = 3 * self.p
        self.vacuum_perm = 4e-7 * np.pi
        self.i_abc = i_abc

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

        # self.total_dofs_bc = len(self.edge_indices)
        self.total_dofs_A_z = len(self.A_z.vector().get_local())
        self.total_dofs_uhat = len(self.uhat.vector().get_local())

        # Partial derivatives in the magnetostatic problem
        self.dR_du = derivative(self.resMS(), self.A_z)
        self.dR_duhat = derivative(self.resMS(), self.uhat)

        # Partial derivatives in the mesh motion subproblem
        self.dRm_duhat = derivative(self.resM(), self.uhat)
        # self.dRm_dedge = self.getBCDerivatives()


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
                self.p,self.s,self.Hc,self.vacuum_perm, self.i_abc)
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

######################## TODO #########################################
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

######################## TODO #########################################
    def setBCMagnetostatic(self):
        """
        Set zero BC on the outer boundary of the domain
        """
        A_outer = Constant(0.0)
        outer_bound = OuterBoundary()
        return DirichletBC(self.V, self.A_z, outer_bound)

    def setBCMeshMotion(self, edge_deltas):
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
            for i in range(self.total_dofs_bc):
                self.uhat_old.vector()[self.edge_indices[i]] += self.increment_deltas[i]
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
        solver_ms.parameters['snes_solver']['linear_solver']='mumps'
        solver_ms.parameters['snes_solver']['error_on_nonconvergence'] = False
        # solver_ms.solve()
        self.B = project(as_vector((self.A_z.dx(1),-self.A_z.dx(0))),
                        VectorFunctionSpace(self.mesh,'DG',0))

        area_func_unit = interpolate(Constant(1.0), self.A_z.function_space())

        subdomain_range = np.arange(41,112+1)

        #### For debugging purposes ###
        for i in subdomain_range:
            print(i, self.getSubdomainArea(area_func_unit, i))
        ###############################

        self.winding_delta_A_z = np.array(
            self.extractSubdomainAverageA_z(
                func=self.A_z,
                subdomain_range=subdomain_range
        )[1])


        self.winding_area = self.getSubdomainArea(
            func_unit=area_func_unit,
            subdomain=42
        )

        self.magnet_area = self.getSubdomainArea(
            func_unit=area_func_unit,
            subdomain=29
        )

        steel_subdomains = [1, 2, 3]
        self.steel_area = 0
        for subdomain in steel_subdomains:
            self.steel_area += self.getSubdomainArea(
                func_unit=area_func_unit,
                subdomain=subdomain
            )
        # NOTES FOR RU:
        # NECESSARY OUTPUTS ARE:
        #   - self.winding_delta_A_z
        #   - self.winding_area
        #   - self.magnet_area
        #   - self.steel_area


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
    iq                  = 282.2 / 3
    i_abc               = [
        -iq * np.sin(0.),
        -iq * np.sin(-2*np.pi/3),
        -iq * np.sin(2*np.pi/3),
    ]
    problem = MotorProblem(i_abc=i_abc)
    # problem.solveMeshMotion()
    problem.solveMagnetostatic()

    vtkfile_A_z = File('solutions/Magnetic_Vector_Potential.pvd')
    vtkfile_B = File('solutions/Magnetic_Flux_Density.pvd')
    vtkfile_A_z << problem.A_z
    vtkfile_B << problem.B

    print("winding area:", problem.winding_area)
    print("magnet area:", problem.magnet_area)
    print("steel area:", problem.steel_area)
    ###### Test the average calculation for the flux linkage
    subdomain_range = np.arange(41,112)
    # asdf, deltas  = problem.extractSubdomainAverageA_z(
    #     func=problem.A_z,
    #     subdomain_range=subdomain_range
    # )
    # for i in range(len(subdomain_range)):
    #     print("Average A_z for subdomain "+str(i+41))
    #     # print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
    #     print(asdf[i])
    #
    # for i in range(len(deltas)):
    #     print("Delta A_z for Stator Tooth "+str(i+1))
    #     # print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
    #     print(deltas[i])

    plt.figure(1)
    plot(problem.A_z)
    # plot(problem.B, linewidth=40)
#    ALE.move(problem.mesh, problem.uhat)
    # plot(problem.mesh)
#    print(problem.A_z.vector().get_local()[:10])
    plt.show()
