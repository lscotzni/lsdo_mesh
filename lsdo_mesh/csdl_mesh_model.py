import numpy as np
import csdl
from csdl import Model

class MeshModel(Model):
    def initialize(self):
        
        self.parameters.declare('ffd_parametrization')
        self.parameters.declare('edge_parametrization')
        self.parameters.declare('mesh_points')
        self.parameters.declare('ffd_cps')
        self.parameters.declare('num_points')

        # STEPS:
        # 1. OPTIMIZER OUTPUTS UPDATED SHAPE PARAMETER VALUES
        # 2. APPLY MAT-VEC PRODUCT BETWEEN SHAPE PARAMETER SPARSE MAT AND SHAPE PARAMETER DELTAS;
        #       - THEN ADD TO ORIGINAL FFD CPs (FOUND IN ffd_cp_instances)
        # 3. APPLY MAT-VEC PRODUCT BETWEEN FFD SPARSE MAT AND NEW FFD CPs
        #       - ADD TO ORIGINAL POINT COORDINATES (FOUND IN mesh_points_instances)
        # 4. APPLY MAT-VEC PRODUCT BETWEEN EDGE PROJECTION AND NEW POINT COORDINATES
        #       - PRODUCES LOCATION OF NEW POINTS

    def define(self):
        
        ffd_parametrization                 = self.parameters['ffd_parametrization']
        edge_parametrization                = self.parameters['edge_parametrization']
        mesh_points                         = self.parameters['mesh_points']
        ffd_cps                             = self.parameters['ffd_cps']
        num_points                          = self.parameters['num_points']

        # SHAPE PARAMETER CONCATENATION:
        # for loop declaring shape parameters as variables
        # use create_output to initialize the empty shape_parameter_vec vector
        # within for loop, assign different shape parameters to parts of empty shape_parameter_vec

        delta_ffd_cp = self.declare_variable('delta_ffd_cp')
        
        # getting deltas of the defined mesh points
        
        delta_mesh_points = csdl.matvec(ffd_parametrization, delta_ffd_cp) 
        
        delta_mesh_points = self.register_output(
            'delta_mesh_points',
            var=delta_mesh_points,
            # shape=delta_mesh_points.shape,
        )

        print(delta_mesh_points)
        print(delta_ffd_cp)
        

            
        mesh_instance = self.declare_variable(
            'mesh_points',
            val=mesh_points_instances[i],
            shape=mesh_points_instances[i].shape
        )
    
        new_mesh_points = mesh_instance + delta_mesh_points
        new_mesh_points = self.register_output(
            'new_mesh_points',
            var=new_mesh_points,
            # shape=new_mesh_points.shape
        )
    
        edge_param_sps_mat = edge_parametrization_instances[i]
        # edge_param_sps_mat = self.declare_variable(
        #     'edge_param_sps_mat_{}'.format(i+1),
        #     val=edge_param_sps_mat,
        #     shape=edge_param_sps_mat.shape
        # )
    

        new_edge_nodes = csdl.matvec(
            edge_param_sps_mat, new_mesh_points
        )

        self.register_output(
            'new_edge_nodes',
            new_edge_nodes,
            # new_edge_nodes.shape
        )

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
        shape_param_vec = []
        for i, name in enumerate(shape_parameter_list_input):
            local_SP = self.declare_variable(name=name) # declaring variable for actual SP
            for j in range(shape_parameter_index_input[i]):
                shape_param_vec.append(local_SP) # appending SP variable "n" times based on 'linear' or 'constant'
                # 'constant': n = 1
                # 'linear': n = 2

        shape_param_vec = self.declare_variable(
            name='shape_param_vec',
            val=shape_param_vec
        )
        
        # shape_param_vec = self.declare_variable(
        #     name='shape_parameter_vec',
        #     shape=(shape_parametrization.shape[1],)
        # )

        delta_ffd_cp = csdl.matvec(shape_parametrization, shape_param_vec)

        delta_ffd_cp = self.register_output(
            'delta_ffd_cp',
            var=delta_ffd_cp,
            # shape=(shape_parametrization.shape[0],),
        )