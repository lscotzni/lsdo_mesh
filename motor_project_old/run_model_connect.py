import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

from motor_fea import *
# MODEL IMPORTS (IN ORDER OF self.add() CALLS IN CSDL)
from motor_mesh_new import motor_mesh_generator

from mesh_motion_model import *
from magnetic_states_model import *

from csdl_fenics_models.flux_influence_ec_loss_model import FluxInfluenceECModel
from csdl_fenics_models.flux_influence_hysteresis_loss_model import FluxInfluenceHModel
from csdl_fenics_models.A_z_air_gap_model import AzAirGapModel
from csdl_fenics_models.area_model import AreaModel

from post_processing_dir.flux_linkage_model import FluxLinkageModel
from post_processing_dir.electrical_model import ElectricalModel

from post_processing_dir.time_average_model import TimeAverageModel

from post_processing_dir.power_loss_model import PowerLossModel
from post_processing_dir.torque_loss_model import TorqueLossModel
from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.mass_model import MassModel


class ConnectionTestModel(Model):
    def initialize(self):
        self.parameters.declare('mesh_objects')
        self.parameters.declare('fea_list')
        self.parameters.declare('electrical_angles')

    def define(self):
        mesh_object = self.parameters['mesh_objects'] # PYTHON OBJECT HOLDING FFD CSDL MODELS
        fea_list    = self.parameters['fea_list']
        elec_angles = self.parameters['electrical_angles']
        instances   = len(fea_list)

        shape_param_model   = vars(mesh_object)['shape_parameter_model']
        mesh_models         = vars(mesh_object)['mesh_models']

        self.add(shape_param_model, 'shape_parameter_model')

        self.declare_variable('num_windings', val=13.)
        
        for i in range(instances):

            instance_model = Model()

            instance_model.declare_variable('num_windings', val=13.)

            instance_model.add(mesh_models[i], 'mesh_model_{}'.format(i+1))
            # NOTE: NEED TO CONNECT edge_deltas FROM shape_parameter_model TO mesh_model_i
            # instance_model.add(MeshMotionModel(fea=fea_list[i]), 'mesh_motion_model_{}'.format(i+1))
            instance_model.add(MagneticStatesModel(fea=fea_list[i]), 'magnetic_states_model_{}'.format(i+1))
            instance_model.add(FluxInfluenceECModel(fea=fea_list[i]), 'flux_influence_ec_model_{}'.format(i+1))
            instance_model.add(FluxInfluenceHModel(fea=fea_list[i]), 'flux_influence_h_model_{}'.format(i+1))
            instance_model.add(AzAirGapModel(fea=fea_list[i]), 'air_gap_A_z_model_{}'.format(i+1))
            instance_model.add(AreaModel(fea=fea_list[i]), 'area_model_{}'.format(i+1))

            instance_model.add(FluxLinkageModel(), 'flux_linkage_model_{}'.format(i+1))
            instance_model.add(ElectricalModel(theta=elec_angles[i]), 'electrical_model_{}'.format(i+1))
            
            self.add(instance_model, name='instance_model_{}'.format(i+1), promotes=[])

        self.add(TimeAverageModel(instances=instances), name='time_average_model')

        for i in range(instances):
            self.connect('instance_model_{}.B_influence_ec'.format(i+1),'B_influence_ec_{}'.format(i+1))
            self.connect('instance_model_{}.B_influence_h'.format(i+1),'B_influence_h_{}'.format(i+1))
            self.connect('instance_model_{}.input_power'.format(i+1),'input_power_{}'.format(i+1))
            self.connect('instance_model_{}.electromagnetic_torque'.format(i+1),'electromagnetic_torque_{}'.format(i+1))
            # self.connect('num_windings', 'instance_model_{}.num_windings'.format(i+1))

        self.add(PowerLossModel(), 'power_loss_model')
        self.add(TorqueLossModel(), 'torque_loss_model')
        self.add(EfficiencyModel(), 'efficiency_model')
        self.add(MassModel(), 'mass_model')

        

        # NOTE-S:
        # ADD AIR GAP MODEL, VECTOR POT MODEL AND ELECTRICAL MODEL TO INSTANCE MODEL
        #   - THE ABOVE 3 MODELS NEED TO BE ALTERED TO FIT THE STYLE (REMOVE INSTANCE IN INDIVIDUAL MODEL, ETC.)
        # TIME AVERAGE MODEL NEEDS TO BE INCLUDED
        # POWER LOSS, TORQUE LOSS, EFFICIENCY AND MASS NEED TO BE INCLUDED



if __name__ == '__main__':
    shift           = 2.5
    mech_angles     = np.arange(0,30,5)
    # rotor_rotations = np.pi/180*np.arange(0,30,5)
    rotor_rotations = np.pi/180*mech_angles[:1]
    instances       = len(rotor_rotations)

    mesh_object     = motor_mesh_generator(
        base_angle=shift*np.pi/180,
        rotation_angles=rotor_rotations,
        file_name='mesh_files/motor_mesh',
        poles=12,
        coarse_test=False,
    )

    fea_list        =  []
    for i in range(instances):
        f = open('edge_deformation_data/init_edge_coords_{}.txt'.format(i+1), 'r+')
        old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
        f.close()

        f = open('edge_deformation_data/edge_coord_deltas_{}.txt'.format(i+1), 'r+')
        edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
        f.close()
        
        f = open('edge_deformation_data/A_z_air_gap_coords_1.txt', 'r+')
        A_z_air_gap_coords = np.fromstring(f.read(), dtype=float, sep=' ')
        f.close()

        fea = MotorFEA(mesh_file='mesh_files/motor_mesh_{}'.format(i+1), old_edge_coords=old_edge_coords)
        fea.edge_deltas         = 0.1 * edge_deltas
        fea.A_z_air_gap_indices = fea.locateAzIndices(A_z_air_gap_coords)

        fea_list.append(fea)
        

    connect_test_model  = ConnectionTestModel(
        mesh_objects=mesh_object,
        fea_list=fea_list,
        electrical_angles=(shift + rotor_rotations) * np.pi/180
    )

    sim     = Simulator(connect_test_model)

    sim.run()

    print(sim['instance_model_1.B_influence_ec'])
    print(sim['instance_model_1.B_influence_h'])