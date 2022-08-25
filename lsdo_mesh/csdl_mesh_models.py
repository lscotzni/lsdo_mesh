import numpy as np
import csdl
from csdl import Model

class EdgeUpdateModel(Model):
    def initialize(self):
        
        self.parameters.declare('ffd_parametrization')
        self.parameters.declare('edge_parametrization')
        self.parameters.declare('mesh_points')
        self.parameters.declare('ffd_cps')
        self.parameters.declare('initial_edge_coords')

        # STEPS:
        # 1. OPTIMIZER OUTPUTS UPDATED SHAPE PARAMETER VALUES
        # 2. APPLY MAT-VEC PRODUCT BETWEEN SHAPE PARAMETER SPARSE MAT AND SHAPE PARAMETER DELTAS;
        #       - THEN ADD TO ORIGINAL FFD CPs (FOUND IN ffd_cp_instances)
        # 3. APPLY MAT-VEC PRODUCT BETWEEN FFD SPARSE MAT AND NEW FFD CPs
        #       - ADD TO ORIGINAL POINT COORDINATES (FOUND IN mesh_points_instances)
        # 4. APPLY MAT-VEC PRODUCT BETWEEN EDGE PROJECTION AND NEW POINT COORDINATES
        #       - PRODUCES LOCATION OF NEW EDGE POINTS

    def define(self):
        
        ffd_parametrization     = self.parameters['ffd_parametrization']
        edge_parametrization    = self.parameters['edge_parametrization']
        mesh_points             = self.parameters['mesh_points']
        ffd_cps                 = self.parameters['ffd_cps']
        initial_edge_coords     = self.parameters['initial_edge_coords']

        # SHAPE PARAMETER CONCATENATION:
        # for loop declaring shape parameters as variables
        # use create_output to initialize the empty shape_parameter_vec vector
        # within for loop, assign different shape parameters to parts of empty shape_parameter_vec

        delta_ffd_cp = self.declare_variable('delta_ffd_cp', shape=(ffd_parametrization.shape[1],))
        # NOTE: SHAPES MUST ALWAYS BE SET UP (WILL NOT AUTOMATICALLY BE READ IN CSDL)
        
        # getting deltas of the defined mesh points
        
        delta_mesh_points = csdl.matvec(ffd_parametrization, delta_ffd_cp) 
        
        delta_mesh_points = self.register_output(
            'delta_mesh_points',
            var=delta_mesh_points,
        )
        
        mesh_instance = self.declare_variable(
            'mesh_points',
            val=mesh_points,
            shape=mesh_points.shape
        )
    
        new_mesh_points = mesh_instance + delta_mesh_points
        new_mesh_points = self.register_output(
            'new_mesh_points',
            var=new_mesh_points,
            # shape=new_mesh_points.shape
        )

        # need edge points here to add to the deltas
        edge_param_sps_mat = edge_parametrization
    

        new_edge_nodes = csdl.matvec(
            edge_param_sps_mat, new_mesh_points
        )

        # edge_deltas_temp = new_edge_nodes[:-2] - initial_edge_coords

        edge_deltas = self.create_output(
            'edge_deltas',
            val=0.,
            shape=initial_edge_coords.shape,
        )

        for i in range(int(initial_edge_coords.shape[0]/2)):
        # Ru: only works for the coarsest mesh
        # for i in range(350,600): (use this for block with the newest CSDL ONLY)
            edge_deltas[2*i] = new_edge_nodes[2*i+1] * csdl.cos(new_edge_nodes[2*i]) - \
                initial_edge_coords[2*i+1] * np.cos(initial_edge_coords[2*i])
            edge_deltas[2*i+1] = new_edge_nodes[2*i+1] * csdl.sin(new_edge_nodes[2*i]) - \
                initial_edge_coords[2*i+1] * np.sin(initial_edge_coords[2*i])
            # need the initial theta and r here to do the deformation step
        
        # ------------------------------------------------------------------------
        
        
        # CSDL OPERATIONS:
        # create_input: need to use this for setting the CSDL objects of the parameters  (immutable)
        # create_output: explicitly computed output (will use this for intermediate calculations)
        # 
        # register_output: we will use this to register the final set of edge node coordinates
        # new_point_location = matvec(ffd_param, ffd_deltas) + original_points

class ShapeParameterModel(Model):
    def initialize(self):
        self.parameters.declare('shape_parameter_list_input')
        self.parameters.declare('shape_parameter_index_input')
        self.parameters.declare('shape_parametrization')

    def define(self):
        shape_parameter_list_input          = self.parameters['shape_parameter_list_input']
        shape_parameter_index_input         = self.parameters['shape_parameter_index_input']
        shape_parametrization               = self.parameters['shape_parametrization']
        # the shape parameter list is coming directly from lsdo_mesh; the user will define the DVs and SPs and connect them,
        # and in here we need to add the SPs as variables and concatenate into a large list looping through the list

        shape_parameter_total   = int(np.dot(np.ones((len(shape_parameter_list_input),)),shape_parameter_index_input))
        shape_param_vec = self.create_output(
            name='shape_param_vec',
            shape=(shape_parameter_total,)
        )

        counter = 0
        for i, name in enumerate(shape_parameter_list_input):
            local_SP    = self.declare_variable(name=name) # declaring variable for actual SP
            len_SP      = shape_parameter_index_input[i]
            # print(len_SP)
            expanded_local_SP  = csdl.expand(local_SP, (len_SP,))
            # print(expanded_local_SP.shape)
            shape_param_vec[counter:counter + len_SP] = expanded_local_SP
            counter += len_SP

        delta_ffd_cp = csdl.matvec(shape_parametrization, shape_param_vec)

        delta_ffd_cp = self.register_output(
            'delta_ffd_cp',
            var=delta_ffd_cp,
        )