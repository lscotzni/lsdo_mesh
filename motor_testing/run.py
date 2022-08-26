import numpy as np
import os
import csdl
from csdl_om import Simulator
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
            - magnet_thickness_sp
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

        magnet_thickness_dv = self.create_input('magnet_thickness_dv', val=0.)
        magnet_thickness_sp = self.register_output(
            'magnet_thickness_sp',
            1*magnet_thickness_dv
        )

        '''
        THE FINAL OUTPUTS HERE ARE THE SHAPE PARAMETERS THAT FEED INTO THE 
        INDIVIDUAL MESH MODELS WITHIN INSTANCE MODELS
        LIST OF SHAPE PARAMETERS:
            - inner_stator_radius_sp
            - magnet_thickness_sp
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
                    mesh_points=param_dict['mesh_points'][i],
                    ffd_cps=param_dict['ffd_cps'][i],
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
    for i in range(num_vertices):
        delta[i, 1] = -.02 # shifting first magnet + air gaps radially in by 0.02 m
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
    mm.create_motor_mesh()
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
    sim['magnet_thickness_dv'] = -0.02
    sim.run()
    # sim.check_partials(compact_print=False)
    # edge_deltas_csdl = sim['edge_deltas']



    ''' --- GETTING DATA FILES FOR INITIAL EDGE COORDS AND EDGE DELTAS BELOW --- '''
    m = mm.motor_mesh_object

    edge_coords_dir = 'edge_deformation_data'
    init_edge_coords    = 'init_edge_coords_'
    edge_coord_deltas   = 'edge_coord_deltas_'

    if coarse_test:
        edge_coords_dir = edge_coords_dir + '_coarse'
        init_edge_coords  = init_edge_coords + 'coarse_'
        edge_coord_deltas = edge_coord_deltas + 'coarse_'

    if not os.path.isdir(edge_coords_dir):
        os.mkdir(edge_coords_dir)

    if os.path.isdir(edge_coords_dir) == False:
        os.mkdir(edge_coords_dir)

    for instance in range(len(rotor_rotations)):
        old_edge_coords = getInitialEdgeCoords(mesh_object=m, instance=instance)
        edge_deltas     = generateMeshMovement(mesh_object=m, instance=instance)

        # init_edge_coords = 'init_edge_coords_{}.txt'.format(i+1)
        # f1  = open(init_edge_coords, 'w')
        # for i in range(old_edge_coords.shape[0]):
        #     f.write(old_edge_coords[i] + '\n')
        # f1.close()
        np.savetxt(
            edge_coords_dir + '/' + init_edge_coords + '{}.txt'.format(instance+1), 
            old_edge_coords
        )

        # edge_coord_deltas = 'edge_coord_deltas_{}.txt'.format(i+1)
        # f2  = open(edge_coord_deltas, 'w')
        # for i in range(edge_deltas.shape[0]):
        #     f.write(edge_deltas[i] + '\n')
        # f2.close()
        np.savetxt(
            edge_coords_dir + '/' + edge_coord_deltas + '{}.txt'.format(instance+1), 
            edge_deltas
        )

        np.savetxt(
            edge_coords_dir + '/' + edge_coord_deltas + 'diff_{}.txt'.format(instance+1), 
            sim['edge_update_model_{}.edge_deltas'.format(instance+1)] - edge_deltas
        )
    
    1