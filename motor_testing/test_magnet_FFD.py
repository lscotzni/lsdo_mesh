import numpy as np
import matplotlib.pyplot as plt
import os
import csdl
# from csdl_om import Simulator
from python_csdl_backend import Simulator
from motor_project.geometry.motor_mesh_class import MotorMesh
from lsdo_mesh.csdl_mesh_models import ShapeParameterModel, EdgeUpdateModel

class ShapeParameterUpdateModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('unique_shape_parameter_list')

    def define(self):
        unique_shape_parameter_list = self.parameters['unique_shape_parameter_list']
        print(unique_shape_parameter_list)
        '''
        COMPUTATION OF MAP BETWEEN DESIGN VARIABLES AND SHAPE PARAMETERS

        LIST OF SHAPE PARAMETERS:
            - inner_stator_radius_sp
            - magnet_pos_delta_sp
            - magnet_width_sp
            - outer_stator_radius_sp
            - rotor_radius_sp
            - shaft_radius_sp
            - stator_tooth_shoe_thickness_sp
            - winding_top_radius_sp
            - winding_width_sp
        '''

        # THE IDEA HERE IS TO REGISTER ALL OF THE SHAPE PARAMETERS WITHIN 
        # unique_shape_parameter_list AS OUTPUTS TO FEED INTO THE FFD MODELS
        # OR WE USE THE SHAPE PARAMETER AS DESIGN VARIABLES 
        # shaft_radius_dv = self.create_input('shaft_radius_dv')
        # shaft_radius_sp = self.register_output(
        #     'shaft_radius_sp',
        #     1*shaft_radius_dv
        # )

        magnet_pos_delta_dv = self.create_input('magnet_pos_delta_dv', val=0.)
        magnet_pos_delta_sp = self.register_output(
            'magnet_pos_delta_sp',
            1*magnet_pos_delta_dv
        )

        # magnet_width_dv = self.create_input('magnet_width_dv', val=0.)
        # magnet_width_sp = self.register_output(
        #     'magnet_width_sp',
        #     1*magnet_width_dv
        # )

        '''
        THE FINAL OUTPUTS HERE ARE THE SHAPE PARAMETERS THAT FEED INTO THE 
        INDIVIDUAL MESH MODELS WITHIN INSTANCE MODELS
        LIST OF SHAPE PARAMETERS:
            - inner_stator_radius_sp
            - magnet_pos_delta_sp
            - magnet_width_sp
            - outer_stator_radius_sp
            - rotor_radius_sp
            - shaft_radius_sp
            - stator_tooth_shoe_thickness_sp
            - winding_top_radius_sp
            - winding_width_sp
        '''

class FFDModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('parametrization_dict')
        self.parameters.declare('instances')

    def define(self):

        param_dict = self.parameters['parametrization_dict']
        unique_sp_list = sorted(set(param_dict['shape_parameter_list_input']))
        instances = self.parameters['instances']

        self.add(
            ShapeParameterUpdateModel(
                unique_shape_parameter_list=unique_sp_list
            ),
            'shape_parameter_update_model'
        )

        self.add(
            ShapeParameterModel(
                shape_parameter_list_input=param_dict['shape_parameter_list_input'],
                shape_parameter_index_input=param_dict['shape_parameter_index_input'],
                shape_parametrization=param_dict['shape_parametrization'],
            ),
            'shape_parameter_model'
        )

        for i in range(instances):
            self.add(
                EdgeUpdateModel(
                    ffd_parametrization=param_dict['ffd_parametrization'][i],
                    edge_parametrization=param_dict['edge_parametrization'][i],
                    initial_edge_coords=param_dict['initial_edge_coordinates'][i],
                ),
                'edge_update_model_{}'.format(i+1),
                promotes=[]
            )
            self.connect('delta_ffd_cp', 'edge_update_model_{}.delta_ffd_cp'.format(i+1))

def getInitialEdgeCoords(mesh_object, instance=0):
    m = mesh_object
    old_edge_coords = m.get_ffd_edge_old_coords(output_type='cartesian', instance=instance)
    # Ru: trim out the origin (x=0,y=0) where there's no nearby (dist<1e-10) nodes in the mesh
    return old_edge_coords[:-2]

def generateMeshMovement(mesh_object, angle=0., instance=0):
    m = mesh_object
    # delta = np.zeros((4 * vars(m)['num_ffd_faces'], 2)) # array of deltas applied to FFD FACES
    # delta[:8, 1] = 0
    # for i in range(4):
    #     delta[2 * i, 0] = angle
    #     delta[2 * i + 1, 0] = -angle
    # #delta[:8, 0] = np.pi/6

    # delta[8:, 1] = 0
    # delta[8:, 0] = 0
    num_vertices = 4 * vars(m)['num_ffd_faces']
    delta = np.zeros((num_vertices, 2))
    # delta = np.zeros((8, 2))
    delta_test_radial = -0.002
    delta_test_azim = -0.02
    for i in range(num_vertices):
        delta[i, 1] = delta_test_radial # shifting first magnet + air gaps radially in by 0.02 m
        if np.mod(i,2) == 0:
            delta[i, 0] = delta_test_azim/2 # shifting first magnet + air gaps radially in by 0.02 m
        else:
            delta[i, 0] = -delta_test_azim/2 # shifting first magnet + air gaps radially in by 0.02 m
    edge_deltas= m.test_ffd_edge_parametrization_polar(delta,   
                                                output_type='cartesian',
                                                instance=instance)
    # print(edge_deltas)
    # exit()
    return edge_deltas[:-2]


if __name__ == '__main__':
    shift           = 2.5
    mech_angles     = np.arange(0,30,5)
    # rotor_rotations = np.pi/180*np.arange(0,30,5)
    rotor_rotations = np.pi/180*mech_angles[:2]
    instances       = len(rotor_rotations)

    coarse_test = True

    if coarse_test:
        mesh_file_dir = 'mesh_files_coarse'
        if not os.path.isdir(mesh_file_dir):
            os.mkdir(mesh_file_dir)
        
        mm = MotorMesh(
            file_name=mesh_file_dir + '/motor_mesh',
            popup=False,
            rotation_angles=rotor_rotations,
            base_angle=shift * np.pi/180,
            test=True
        )
    else:
        mesh_file_dir = 'mesh_files'
        if not os.path.isdir(mesh_file_dir):
            os.mkdir(mesh_file_dir)

        mm = MotorMesh(
            file_name=mesh_file_dir + '/motor_mesh',
            popup=False,
            rotation_angles=rotor_rotations,
            base_angle=shift * np.pi/180,
        )
    
    mm.baseline_geometry=True
    mm.magnet_shift_only = True
    mm.create_motor_mesh()
    m = mm.motor_mesh_object
    parametrization_dict    = mm.ffd_param_dict # dictionary holding parametrization parameters
    unique_sp_list = sorted(set(parametrization_dict['shape_parameter_list_input']))
    print(unique_sp_list)

    ''' FFD MODEL '''
    print('Starting FFD CSDL Model')
    ffd_connection_model = FFDModel(
        parametrization_dict=parametrization_dict,
        instances=instances
    )
    # rep = csdl.GraphRepresentation(ffd_connection_model)
    # sim = Simulator(rep)
    sim = Simulator(ffd_connection_model)
    sim['magnet_pos_delta_dv'] = -0.002
    # sim['magnet_width_dv'] = -0.05
    sim.run()
    # sim.visualize_implementation()
    # sim.check_partials(compact_print=False)
    # edge_deltas_csdl = sim['edge_deltas']

    ''' --- GETTING DATA FILES FOR INITIAL EDGE COORDS AND EDGE DELTAS BELOW --- '''
    

    edge_coords_dir = 'edge_deformation_data'
    init_edge_coords    = 'init_edge_coords_'
    edge_coord_deltas   = 'edge_coord_deltas_'

    if coarse_test:
        edge_coords_dir = edge_coords_dir + '_coarse'
        init_edge_coords  = init_edge_coords + 'coarse_'
        edge_coord_deltas = edge_coord_deltas + 'coarse_'

    # if not os.path.isdir(edge_coords_dir):
    #     os.mkdir(edge_coords_dir)

    if os.path.isdir(edge_coords_dir) == False:
        os.mkdir(edge_coords_dir)

    for instance in range(len(rotor_rotations)):
        old_edge_coords = getInitialEdgeCoords(mesh_object=m, instance=instance)
        edge_deltas     = generateMeshMovement(mesh_object=m, instance=instance)

        # SAVING OLD EDGE COORDINATES TO FILE
        np.savetxt(
            edge_coords_dir + '/' + init_edge_coords + '{}.txt'.format(instance+1), 
            old_edge_coords
        )
        # SAVING EDGE DELTAS TO FILE
        np.savetxt(
            edge_coords_dir + '/' + edge_coord_deltas + '{}.txt'.format(instance+1), 
            edge_deltas
        )
        # SAVING DIFFERENCE IN EDGE DELTAS TO FILE
        edge_deltas_csdl = sim['edge_update_model_{}.edge_deltas'.format(instance+1)]
        np.savetxt(
            edge_coords_dir + '/' + edge_coord_deltas + 'diff_{}.txt'.format(instance+1), 
            edge_deltas_csdl - edge_deltas
        )
        # PLOTTING OF OLD EDGE COORDINATES AND EDGE DELTAS FROM METHOD ABOVE
        num_points = int(len(old_edge_coords)/2)
        old_edge_coords_reformat = np.zeros((num_points,2))
        edge_deltas_reformat = np.zeros_like(old_edge_coords_reformat)
        for i in range(num_points):
            old_edge_coords_reformat[i,:] = [old_edge_coords[2*i], old_edge_coords[2*i+1]]
            edge_deltas_reformat[i,:] = [edge_deltas[2*i], edge_deltas[2*i+1]]

        new_edge_coords_reformat = old_edge_coords_reformat + edge_deltas_reformat

        # plt.figure(3*instance + 1)
        # plt.plot(new_edge_coords_reformat[:,0], new_edge_coords_reformat[:,1], 'r*', label='new edge coords')
        # plt.plot(old_edge_coords_reformat[:,0], old_edge_coords_reformat[:,1], 'k*', label='old edge coords')
        # plt.axis('equal')
        # plt.title('Edge deformations from sample methods for instance {}'.format(instance+1))
        # plt.grid()
        # plt.legend()

        # PLOTTING OF OLD EDGE COORDINATES AND EDGE DELTAS CALCULATED FROM CSDL
        num_csdl_points = int(len(edge_deltas_csdl)/2)
        old_edge_coords_reformat_csdl = np.zeros((num_csdl_points,2))
        edge_deltas_reformat_csdl = np.zeros_like(old_edge_coords_reformat_csdl)
        old_edge_coords_csdl = parametrization_dict['initial_edge_coordinates'][instance]
        for i in range(num_csdl_points):
            old_edge_coords_reformat_csdl[i,0] = old_edge_coords_csdl[2*i+1]*np.cos(old_edge_coords_csdl[2*i])
            old_edge_coords_reformat_csdl[i,1] = old_edge_coords_csdl[2*i+1]*np.sin(old_edge_coords_csdl[2*i])
            edge_deltas_reformat_csdl[i,:] = [edge_deltas_csdl[2*i], edge_deltas_csdl[2*i+1]]

        new_edge_coords_reformat_csdl = old_edge_coords_reformat_csdl + edge_deltas_reformat_csdl

        # plt.figure(3*instance + 2)
        # plt.plot(new_edge_coords_reformat_csdl[:,0], new_edge_coords_reformat_csdl[:,1], 'r*', label='csdl new edge coords')
        # plt.plot(old_edge_coords_reformat_csdl[:,0], old_edge_coords_reformat_csdl[:,1], 'k*', label='csdl old edge coords')
        # plt.axis('equal')
        # plt.title('Edge deformations from CSDL FFD Models for instance {}'.format(instance+1))
        # plt.grid()
        # plt.legend()

        plt.figure(3*instance + 3)
        for i in range(new_edge_coords_reformat_csdl.shape[0]):
            if any(edge_deltas_reformat_csdl[i,:]):
                plt.plot(
                    [old_edge_coords_reformat_csdl[i,0], new_edge_coords_reformat_csdl[i,0]],
                    [old_edge_coords_reformat_csdl[i,1], new_edge_coords_reformat_csdl[i,1]],
                    'r-*'
                )
                plt.plot(
                    [0.0, new_edge_coords_reformat_csdl[i,0]],
                    [0.0, new_edge_coords_reformat_csdl[i,1]],
                    'r-'
                )
        plt.plot(old_edge_coords_reformat_csdl[:,0], old_edge_coords_reformat_csdl[:,1], 'k*', label='csdl old edge coords')
        plt.axis('equal')
        plt.title('Edge deformations from CSDL FFD Models for instance {}'.format(instance+1))
        plt.grid()
        plt.legend()

    plt.show()
        