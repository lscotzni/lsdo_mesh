"""
The FEniCS wrapper for the variational forms and the partial derivatives computation
in the magnetostatic problem for the motor on a deformable mesh.
"""
from fea_utils import *
from motor_problem import *
import numpy as np
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh
from scipy.sparse import csr_matrix


class MotorFEA(object):
    """
    The class of the FEniCS wrapper for the motor problem,
    with methods to compute the variational forms, partial derivatives,
    and solve the nonlinear/linear subproblems.
    """
    def __init__(self, mesh_file="motor_mesh_1", i_abc=[0., 0., 0.], old_edge_coords=None):

        self.mesh_file = mesh_file
        # Import the initial mesh from the mesh file to FEniCS
        self.initMesh()
        # Define the function spaces based on the initial mesh
        self.initFunctionSpace()
        # Get the indices of the vertices that would move during optimization
        self.edge_indices = self.locateMeshMotion(old_edge_coords)

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

        self.total_dofs_bc = len(self.edge_indices)
        self.total_dofs_A_z = len(self.A_z.vector().get_local())
        self.total_dofs_uhat = len(self.uhat.vector().get_local())

        # Partial derivatives in the magnetostatic problem
        self.dR_du = derivative(self.resMS(), self.A_z)
        self.dR_duhat = derivative(self.resMS(), self.uhat)

        # Partial derivatives in the mesh motion subproblem
        self.dRm_duhat = derivative(self.resM(), self.uhat)
        


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
        self.winding_id = [42,]
        self.magnet_id = [29,]
        self.steel_id = [1,2,3]
        self.winding_range = range(41,112+1)


    def initFunctionSpace(self):
        """
        Preprocessor 2 to define the function spaces for the mesh motion (VHAT)
        and the problem solution (V)
        """
        self.V = FunctionSpace(self.mesh, 'P', 1)
        self.VHAT = VectorFunctionSpace(self.mesh, 'P', 1)

    def locateMeshMotion(self,old_edge_coords):
        """
        Find the indices of the dofs for setting up the boundary condition
        in the mesh motion subproblem
        """
        V0 = FunctionSpace(self.mesh, 'CG', 1)
        coordinates = V0.tabulate_dof_coordinates()

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

    def areaForm(self, subdomain):
        """
        The UFL form for area calculation of a subdomain
        """
        return Constant(1.0)*J(self.uhat)*self.dx(subdomain)
    
    def funcIntegralForm(self, func, subdomain):
        """
        The UFL form for function integral calculation on a subdomain
        """   
        func_unit = interpolate(Constant(1.0), func.function_space())
        return inner(func, func_unit)*J(self.uhat)*self.dx(subdomain)
        
    def getSubdomainArea(self, subdomains):
        """
        Compute the subdomain area based on its flag
        """
        if type(subdomains) == int:
            subdomain_group = [subdomains]
        else:
            subdomain_group = subdomains
        area = 0
        for subdomain_id in subdomain_group:
            area += self.areaForm(subdomain=subdomain_id)
        return area
        
    def getFuncAverageSubdomain(self, func, subdomain):
        """
        Compute the average function value over a subdomain
        """
        integral = assemble(self.funcIntegralForm(func, subdomain))
        area = assemble(self.getSubdomainArea(subdomain))
        avg_func = integral/area
        return avg_func

    def extractSubdomainAverageA_z(self, func, subdomain_range):
        """
        Get averaged function values over a range of subdomains and 
        the differences between the average values of two neighboring
        subdomains
        """
        subdomain_avg_A_z = []
        subdomain_avg_A_z_deltas = []
        # winding_start, winding_end = 41, 112
        for subdomain in subdomain_range:
            subdomain_avg_A_z.append(
                self.getFuncAverageSubdomain(func, subdomain)
            )
        for i in range(int(len(subdomain_range)/2)):
            subdomain_avg_A_z_deltas.append(
                abs(subdomain_avg_A_z[2*i] - subdomain_avg_A_z[2*i+1])
            )
        return subdomain_avg_A_z, subdomain_avg_A_z_deltas
    
    def calcAreas(self):
        self.winding_area = self.getSubdomainArea(self.winding_id)
        self.magnet_area = self.getSubdomainArea(self.magnet_id)
        self.steel_area = self.getSubdomainArea(self.steel_id)

        return (assemble(self.winding_area), 
                assemble(self.magnet_area), 
                assemble(self.steel_area))
               

        
    def calcAreaDerivatives(self):
        dWAduhat = assemble(derivative(self.winding_area, self.uhat))
        dMAduhat = assemble(derivative(self.magnet_area, self.uhat))
        dSAduhat = assemble(derivative(self.steel_area, self.uhat))
        
        return (dWAduhat.get_local(), 
                dMAduhat.get_local(), 
                dSAduhat.get_local())
    
            
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
        if min_STEPS >= STEPS:
            STEPS = min_STEPS
        increment_deltas = edge_deltas/STEPS
        return STEPS, increment_deltas

    def setUpMeshMotionSolver(self, report):
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

        # Nonlinear solver parameters
        REL_TOL_M = 1e-4
        MAX_ITERS_M = 10

        # Set up nonlinear problem for mesh motion:
        problem_m = NonlinearVariationalProblem(res_m, self.uhat,
                                                bc_m, Dres_m)
        self.solver_m = NonlinearVariationalSolver(problem_m)
        self.solver_m.parameters['nonlinear_solver']='snes'
        self.solver_m.parameters['snes_solver']['line_search'] = 'bt'
        self.solver_m.parameters['snes_solver']['relative_tolerance'] = REL_TOL_M
        self.solver_m.parameters['snes_solver']['maximum_iterations'] = MAX_ITERS_M
        self.solver_m.parameters['snes_solver']['linear_solver']='mumps'
        self.solver_m.parameters['snes_solver']['error_on_nonconvergence'] = False
        self.solver_m.parameters['snes_solver']['report'] = report

    def solveMeshMotion(self, report=False):

        """
        Solve for the mesh motion subproblem
        """
        # Set up the Newton solver based on the boundary conditions
        self.setUpMeshMotionSolver(report)

        # Incrementally set the BCs to increase to `edge_deltas`
        print(80*"=")
        print(' FEA: total steps for mesh motion:', self.STEPS)
        print(80*"=")
        for i in range(self.STEPS):
            if report == True:
                print(80*"=")
                print("  FEA: Step "+str(i+1)+" of mesh movement")
                print(80*"=")
            self.uhat_old.assign(self.uhat)

            for i in range(self.total_dofs_bc):
                self.uhat_old.vector()[self.edge_indices[i]] += self.increment_deltas[i]
            self.solver_m.solve()
        
        if report == True:
            print(80*"=")
            print(' FEA: L2 error of the mesh motion on the edges:',
                        np.linalg.norm(self.uhat.vector()[self.edge_indices]
                                         - self.edge_deltas))
            print(80*"=")

    def solveMagnetostatic(self, report=False):
        """
        Solve the magnetostatic problem with the mesh movements `uhat`
        """
        if report == True:
            print(80*"=")
            print(" FEA: Solving the magnetostatic problem on the deformed mesh")
            print(80*"=")
        bc_ms = self.setBCMagnetostatic()
        res_ms = self.resMS()
        A_z_ = TrialFunction(self.V)
        Dres_ms = derivative(res_ms, self.A_z, A_z_)

        # Nonlinear solver parameters
        ABS_TOL_M = 1e-6
        REL_TOL_M = 1e-6
        MAX_ITERS_M = 100

        problem_ms = NonlinearVariationalProblem(res_ms, self.A_z,
                                                bc_ms, Dres_ms)
        solver_ms = NonlinearVariationalSolver(problem_ms)
        solver_ms.parameters['nonlinear_solver']='snes'
        solver_ms.parameters['snes_solver']['line_search'] = 'bt'
        solver_ms.parameters['snes_solver']['absolute_tolerance'] = ABS_TOL_M
        solver_ms.parameters['snes_solver']['relative_tolerance'] = REL_TOL_M
        solver_ms.parameters['snes_solver']['maximum_iterations'] = MAX_ITERS_M
        solver_ms.parameters['snes_solver']['linear_solver']='mumps'
        solver_ms.parameters['snes_solver']['error_on_nonconvergence'] = False
        solver_ms.parameters['snes_solver']['report'] = report
        solver_ms.solve()
        self.B = project(as_vector((self.A_z.dx(1),-self.A_z.dx(0))),
                        VectorFunctionSpace(self.mesh,'DG',0))

        self.winding_delta_A_z = np.array(
         self.extractSubdomainAverageA_z(
             func=self.A_z,
             subdomain_range=self.winding_range
        )[1])

    # TODO
    def calcWindingAz(self):
        winding_A_z = []
        for subdomain in self.winding_range:
            winding_A_z.append(
                self.getFuncAverageSubdomain(self.A_z, subdomain)
            )
        return winding_A_z
    
    # TODO
    def calcWindingAzDerivatives(self):
        dAz = []
        duhat = []
        for subdomain in self.winding_range:
            func_integral_form = self.funcIntegralForm(self.A_z, subdomain)
            area_form = self.areaForm(subdomain)
            func_integral_i = assemble(func_integral_form)
            area_i = assemble(area_form)
            dAz_i = 1/area_i*assemble(derivative(func_integral_form, self.A_z))
            
            # quotient rule for derivatives; Nu: numerator, De: denominator
            dNudu = assemble(derivative(func_integral_form, self.uhat))
            dDedu = assemble(derivative(area_form, self.uhat))
            
            duhat_i = (area_i*dNudu + func_integral_i*dDedu)/area_i**2
            
            
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
    f = open('edge_deformation_data/init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('edge_deformation_data/edge_coord_deltas.txt', 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

#    print("Number of nonzero displacements:", np.count_nonzero(edge_deltas))
    
    # One-time computation for the initial edge coordinates from
    # the code that creates the mesh file
    # old_edge_coords = getInitialEdgeCoords()
    
    problem = MotorFEA(mesh_file="mesh_files/motor_mesh_1", i_abc=i_abc, 
                            old_edge_coords=old_edge_coords)
    
    problem.edge_deltas = edge_deltas
    problem.solveMeshMotion(report=True)
#    plt.figure(1)
#    problem.moveMesh(problem.uhat)
#    plot(problem.mesh)
#    plt.show()

    problem.solveMagnetostatic(report=True)
#    problem.calcWindingAzDerivatives()
#    vtkfile_A_z = File('solutions/Magnetic_Vector_Potential.pvd')
#    vtkfile_B = File('solutions/Magnetic_Flux_Density.pvd')
#    vtkfile_A_z << problem.A_z
#    vtkfile_B << problem.B

#    
#    ###### Test the average calculation for the flux linkage
#    subdomain_range = range(41,112)
#    asdf, deltas  = problem.extractSubdomainAverageA_z(
#        func=problem.A_z,
#        subdomain_range=subdomain_range
#    )
#    
#    for i in range(len(subdomain_range)):
#        print("Average A_z for subdomain "+str(i+41))
##        print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
#        print(asdf[i])

#    for i in range(len(deltas)):
#        print("Delta A_z for Stator Tooth "+str(i+1))
#        # print(problem.getFuncAverageSubdomain(func=problem.A_z, subdomain=i+1))
#        print(deltas[i])
#    print('------')

##    for i, ind in enumerate(subdomain_range):
##        func=problem.A_z
##        func_unit = interpolate(Constant(1.0), func.function_space())
##        print("Area of Stator Winding "+str(ind))
##        print(problem.getSubdomainArea(subdomain=ind))

#    # for i, ind in enumerate(subdomain_range):
#    #     func=problem.A_z
#    #     func_unit = interpolate(Constant(1.0), func.function_space())
#    #     print("Area of Stator Winding"+str(ind))
#    #     print(problem.getSubdomainArea(subdomain=i+1))
##    
#    plt.figure(1)
#    asdf = plot(problem.A_z)
#    plt.colorbar(asdf)
#    # plot(problem.B, linewidth=40)
#    # ALE.move(problem.mesh, problem.uhat)
#    # plot(problem.mesh)
#    # print(problem.A_z.vector().get_local()[:10])
#    plt.show()
