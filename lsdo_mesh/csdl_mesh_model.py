import numpy as np
import csdl
from csdl import Model

class MeshModel(Model):
    def initialize(self):
        self.parameters.declare('shape_parametrization')
        self.parameters.declare('ffd_parametrization')
        self.parameters.declare('edge_parametrization_instances')
        self.parameters.declare('mesh_points_instances')
        self.parameters.declare('ffd_cp_instances')
        self.parameters.declare('num_mesh_instances')
        self.parameters.declare('num_points')

        # made initially but not sure if this is the correct method
        # shape_sps_mat = self.create_input(
        #     name='shape_sps_mat',
        #     val=shape_parametrization,
        #     shape=(shape_parametrization.shape[0], shape_parametrization.shape[1])
        # )

        # ffd_sps_mat = self.create_input(
        #     name='ffd_sps_mat',
        #     val=ffd_parametrization,
        #     shape=ffd_parametrization.shape
        # )

        # edge_sps_mat = self.create_input(
        #     name='edge_sps_mat',
        #     val=edge_parametrization,
        #     shape=edge_parametrization.shape
        # )

        # mesh_points_instances = self.create_input(
        #     name='mesh_points_instances',
        #     val=mesh_points_instances,
        #     shape=(len(mesh_points_instances),)
        # ) # size (n,1) for n instances, where each entry holds a numpy array
        
        # ffd_control_points = self.create_input(
        #     name='ffd_control_points_instances',
        #     val=ffd_cp_instances,
        #     shape=(len(ffd_cp_instances),)
        # ) # size (n,1) for n instances, where each entry holds a numpy array


    def define(self):
        shape_parametrization               = self.parameters['shape_parametrization']
        ffd_parametrization                 = self.parameters['ffd_parametrization']
        edge_parametrization_instances      = self.parameters['edge_parametrization_instances']
        mesh_points_instances               = self.parameters['mesh_points_instances']
        ffd_cp_instances                    = self.parameters['ffd_cp_instances']
        num_mesh_instances                  = self.parameters['num_mesh_instances']
        num_points                          = self.parameters['num_points']

        # print(mesh_points_instances)
        # print(edge_parametrization_instances)
        # print(len(mesh_points_instances))
        # print(len(mesh_points_instances[0]))
        # print(len(mesh_points_instances[1]))
        # print(len(edge_parametrization_instances))
        # print(edge_parametrization_instances[0].shape)
        # print(edge_parametrization_instances[1].shape)
        # exit()
        
        # edge_sps_mat = self.create_input(
        #     name='edge_sps_mat',
        #     val=edge_parametrization,
        #     shape=edge_parametrization.shape
        # )

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
        
        for i in range(mesh_instances):
            
            mesh_instance = self.declare_variable(
                'mesh_instance_{}'.format(i+1),
                val=mesh_points_instances[i],
                shape=mesh_points_instances[i].shape
            )
        
            new_mesh_points = mesh_instance + delta_mesh_points
            print(type(new_mesh_points))
            new_mesh_points = self.register_output(
                'new_mesh_points_{}'.format(i+1),
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
                'new_edge_nodes_{}'.format(i+1),
                new_edge_nodes,
                # new_edge_nodes.shape
            )
                # print(new_edge_nodes.shape)
            # exit()
        # ------------------------------------------------------------------------
        
        # STEPS:
        # 1. OPTIMIZER OUTPUTS UPDATED SHAPE PARAMETER VALUES
        # 2. APPLY MAT-VEC PRODUCT BETWEEN SHAPE PARAMETER SPARSE MAT AND SHAPE PARAMETER DELTAS;
        #       - THEN ADD TO ORIGINAL FFD CPs (FOUND IN ffd_cp_instances)
        # 3. APPLY MAT-VEC PRODUCT BETWEEN FFD SPARSE MAT AND NEW FFD CPs
        #       - ADD TO ORIGINAL POINT COORDINATES (FOUND IN mesh_points_instances)
        # 4. APPLY MAT-VEC PRODUCT BETWEEN EDGE PROJECTION AND NEW POINT COORDINATES
        #       - PRODUCES LOCATION OF NEW POINTS
        # CSDL OPERATIONS:
        # create_input: need to use this for setting the CSDL objects of the parameters  (immutable)
        # create_output: explicitly computed output (will use this for intermediate calculations)
        # 
        # register_output: we will use this to register the final set of edge node coordinates
        # new_point_location = matvec(ffd_param, ffd_deltas) + original_points