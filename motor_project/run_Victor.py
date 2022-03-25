import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

from motor_mesh_new import MotorMeshGenerator
from post_processing import PostProcessingModel
from motor_fea import *

# -------------------- MOTOR MODEL --------------------
class MotorModel(Model):
    def initialize(self):
        self.parameters.declare('shape_param_object')
        self.parameters.declare('mesh_objects')
        self.parameters.declare('motor_length')
        self.parameters.declare('omega')
        self.parameters.declare('current_amplitude')
        self.parameters.declare('fea_list')
        self.parameters.declare('frequency')
        self.parameters.declare('angles')

    def define(self):
        shape_param_model   = self.parameters['shape_param_object']
        mesh_objects        = self.parameters['mesh_objects']
        motor_length        = self.parameters['motor_length']
        omega               = self.parameters['omega']
        frequency           = self.parameters['frequency']
        current_amplitude_i = self.parameters['current_amplitude']
        fea_list            = self.parameters['fea_list']
        angles              = self.parameters['angles']
        instances           = len(angles)

        self.add(shape_param_model, name='shape_param_model')

        self.declare_variable('delta_ffd_cp')

        for i in range(instances):
            self.add(
                mesh_objects[i],
                name='mesh_motion_model_{}'.format(i+1)
            )

        self.add(
            PostProcessingModel(
                motor_length=motor_length,
                omega=omega,
                current_amplitude=iq,
                fea_list=fea_list,
                frequency=freq,
                angles=np.array(angles[:instances]) * p / 2
            ),
            name='post_processing_model'
        )

# -------------------- CODE --------------------
shift = 2.5
 
# rotor_rotations = np.pi*180*np.arange(0,30,5)
rotor_rotations = np.pi*180*np.arange(0,5,5)

mesh_objects = MotorMeshGenerator(
    base_angle=np.pi / 180 * shift,
    # rotation_angles=np.pi / 180 * np.array([0]), 
    rotation_angles=rotor_rotations, 
    file_name='mesh_files/motor_mesh',
    poles=12,
    coarse_test=False
)  
print(vars(mesh_objects)['mesh_model'])
print(vars(mesh_objects)['shape_parameter_model'])

slot_fill_factor    = 0.6
num_windings        = 13
rpm                 = 1000.
iq                  = 20.

p                   = 12
t                   = 0
motor_length        = 70.e-3
angle_shift         = 2.5
angles              = (angle_shift + np.arange(0,30,5)) * np.pi / 180
# instances           = len(angles)
instances           = len(rotor_rotations)
omega               = rpm * 2 * np.pi / 60.
freq                = rpm * p / 120.

# RUNNING FEA FOR EACH ROTATION INSTANCE
fea_list = []
for i in range(instances):
    print('-------')
    print('fea instance:', i)
    print('-------')
    # angle = 0.0 + i * 5 * np.pi / 180
    angle = angles[i] * p / 2
    # i_abc               = [
    #     iq * np.sin(angle),
    #     iq * np.sin(angle + 2*np.pi/3),
    #     iq * np.sin(angle - 2*np.pi/3),
    # ]
    f = open('edge_deformation_data/init_edge_coords_{}.txt'.format(i+1), 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('edge_deformation_data/edge_coord_deltas_{}.txt'.format(i+1), 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    fea = MotorFEA(mesh_file="mesh_files/motor_mesh_{}".format(i+1),
                        old_edge_coords=old_edge_coords)
    fea.angle = angle
    updateR(fea.iq, iq / (0.00016231 / num_windings))
    fea.solveMagnetostatic(report=True)
    

    f = open('edge_deformation_data/A_z_air_gap_coords_{}.txt'.format(i+1), 'r+')
    A_z_air_gap_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()
    
    fea.edge_deltas = 0.1*edge_deltas # doesn't matter right now
    fea.A_z_air_gap_indices = fea.locateAzIndices(A_z_air_gap_coords)
    print('A_z air gap values: ', fea.extractAzAirGap())
    fea_list.append(fea)

    # SAVING DATA FOR TESTING
    if False:
        vtkfile_A_z = File('post_proc_viz/Magnetic_Vector_Potential_20A_{}.pvd'.format(i))
        vtkfile_B = File('post_proc_viz/Magnetic_Flux_Density_20A_{}.pvd'.format(i))
        vtkfile_A_z << fea.A_z
        vtkfile_B << fea.B

# for i in range(instances):
#     sim['A_z_air_gap_model_{}.A_z'.format(i)] = fea_list[i].A_z.vector().get_local()
#     sim['flux_influence_ec_model_{}.A_z'.format(i)] = fea_list[i].A_z.vector().get_local()
#     sim['flux_influence_ec_model_{}.uhat'.format(i)] = fea_list[i].uhat.vector().get_local()
#     sim['flux_influence_h_model_{}.A_z'.format(i)] = fea_list[i].A_z.vector().get_local()
#     sim['flux_influence_h_model_{}.uhat'.format(i)] = fea_list[i].uhat.vector().get_local()


motor_model     = MotorModel(
    shape_param_object=vars(mesh_objects)['shape_parameter_model'],
    mesh_objects=vars(mesh_objects)['mesh_model'],
    motor_length=motor_length,
    omega=omega,
    current_amplitude=iq,
    fea_list=fea_list,
    frequency=freq,
    angles=np.array(angles[:instances]) * p / 2
)

# motor_model.visualize_sparsity()

# aaa     = sim(motor_model)

# --------------------------- VICTOR'S CSDL DATA ---------------------------
from csdl.utils.graph import min_in_degree, min_out_degree, max_in_degree, max_out_degree, is_tree, count_combined_operations, count_constraints, count_design_variables, count_operations, count_outputs, count_std_operations, count_variables

# set up kwargs for constructing your model if necessary

for x in [0, 1]:
    # create instance of model
    # remember to include necessary kwargs
    m = motor_model

    if x == 0:
        m.optimize_ir(False)
        print('not optimized')
    else:
        print('')
        print('optimized')

    # run compiler front end and middle end manually
    m.define()

    # graph properties
    print('min_in_degree', min_in_degree(m))
    print('min_out_degree', min_out_degree(m))
    print('max_in_degree', max_in_degree(m))
    print('max_out_degree', max_out_degree(m))
    print('is_tree', is_tree(m))

    # operations
    print('num. operations', count_operations(m))
    print('num. std operations', count_std_operations(m))
    print('num. elementwise operations',
          count_std_operations(m, elementwise_only=True))
    print('num. combined operations', count_combined_operations(m))

    # variables
    print('num. variables', count_variables(m))
    print('num. design_variables', count_design_variables(m))
    print('num. constraints', count_constraints(m))
    print('num. outputs', count_outputs(m))

    print('num. variables (not vectorized)',
          count_variables(m, vectorized=False))
    print('num. design_variables (not vectorized)',
          count_design_variables(m, vectorized=False))
    print('num. constraints (not vectorized)',
          count_constraints(m, vectorized=False))
    print('num. outputs (not vectorized)', count_outputs(m, vectorized=False))