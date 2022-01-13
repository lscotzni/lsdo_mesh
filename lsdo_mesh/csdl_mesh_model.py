import numpy as np
import csdl
from csdl import Model

class MeshModel(Model):
    def initialize(self):
        self.parameters.declare('shape_parametrization')
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
        shape_parametrization               = self.parameters['shape_parametrization']
        ffd_parametrization                 = self.parameters['ffd_parametrization']
        edge_parametrization                = self.parameters['edge_parametrization']
        mesh_points                         = self.parameters['mesh_points']
        ffd_cps                             = self.parameters['ffd_cps']
        num_points                          = self.parameters['num_points']

        # SHAPE PARAMETER CONCATENATION:
        # for loop declaring shape parameters as variables
        # use create_output to initialize the empty shape_parameter_vec vector
        # within for loop, assign different shape parameters to parts of empty shape_parameter_vec

        shape_param_vec = self.declare_variable(
            name='shape_parameter_vec',
            shape=(shape_parametrization.shape[1],)
        )

        # getting deltas of the FFD control points
        delta_ffd_cp = csdl.matvec(shape_parametrization, shape_param_vec)
        # FORMAT TEST BELOW
        # delta_ffd_cp = self.register_output(
        #     'delta_ffd_cp',
        #     var=delta_ffd_cp,
        #     # shape=(shape_parametrization.shape[0],),
        # )

        # NOT NEEDED; 
        # delta_ffd_cp = self.declare_variable(
        #     'delta_ffd_cp',
        #     val=delta_ffd_cp,
        #     shape=(shape_parametrization.shape[0],),
        # )
        
        # getting deltas of the defined mesh points
        print(delta_ffd_cp)
        delta_mesh_points = csdl.matvec(ffd_parametrization, delta_ffd_cp) 
        print(delta_mesh_points)
        # delta_mesh_points = self.register_output(
        #     'delta_mesh_points',
        #     var=delta_mesh_points,
        #     # shape=delta_mesh_points.shape,
        # )
        

            
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