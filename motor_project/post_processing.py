import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
import time

from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.flux_linkage_model import FluxLinkageModel
from post_processing_dir.electrical_model import ElectricalModel
from post_processing_dir.power_loss_model import PowerLossModel
from post_processing_dir.mass_model import MassModel
from post_processing_dir.vector_potential_model import VectorPotentialModel


from B_rms_2_model import Brms2Model
from A_z_air_gap_model import AzAirGapModel
from area_model import AreaModel

from motor_fea import *

class PostProcessingModelDuplicate(Model):
    def define(self):
        flux_linkage_model  = self.add(FluxLinkageModel(), name='flux_linkage_model')
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        # final form should look like this
        # variables like winding_area will be defined as CSDL variables in the FEA
        # file, so they don't need to be declared here

class PostProcessingModel(Model):
    def initialize(self):
        self.parameters.declare('winding_delta_A_z')
        self.parameters.declare('winding_area')
        self.parameters.declare('motor_length')
        self.parameters.declare('omega')
        self.parameters.declare('magnet_area')
        self.parameters.declare('steel_area')
        self.parameters.declare('phase_current_dq')
        self.parameters.declare('fea')
        self.parameters.declare('frequency')

    def define(self):
        # winding_delta_A_z   = self.parameters['winding_delta_A_z']
        # winding_area        = self.parameters['winding_area']
        # magnet_area         = self.parameters['magnet_area']
        # steel_area          = self.parameters['steel_area']

        motor_length        = self.parameters['motor_length']
        omega               = self.parameters['omega']
        frequency           = self.parameters['frequency']
        phase_current_dq    = self.parameters['phase_current_dq']
        fea                 = self.parameters['fea']

        # winding_delta_A_z   = self.declare_variable(
        #     name='winding_delta_A_z',
        #     val=winding_delta_A_z,
        #     shape=(36,)
        # )
        # winding_area        = self.declare_variable(
        #     name='winding_area',
        #     val=winding_area
        # )
        # magnet_area         = self.declare_variable(
        #     name='magnet_area',
        #     val=magnet_area
        # )
        # steel_area          = self.declare_variable(
        #     name='steel_area',
        #     val=steel_area
        # )

        motor_length        = self.declare_variable(
            name='motor_length',
            val=motor_length
        )
        omega               = self.declare_variable(
            name='omega',
            val=omega
        )

        frequency           = self.declare_variable(
            name='frequency',
            val=frequency
        )

        phase_current_dq    = self.declare_variable(
            name='phase_current_dq',
            val=phase_current_dq,
            shape=(2,),
        )

        current_amplitude   = self.register_output(
            name='current_amplitude',
            var=phase_current_dq[1]
        )

        num_windings = self.declare_variable(name='num_windings', shape=(1,), val=13.)
        
        # RU'S FEniCS MODELS

        area_model          = self.add(AreaModel(fea=fea), name='area_model')
        B_rms_2_model       = self.add(Brms2Model(fea=fea), name='B_rms_2_model')
        A_z_air_gap_model   = self.add(AzAirGapModel(fea=fea), name='A_z_air_gap_model')

        # SIMPLE P-P MODELS

        vector_pot_model    = self.add(VectorPotentialModel(), name='vector_pot_model')
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        power_loss_model    = self.add(PowerLossModel(), name='power_loss_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        mass_model          = self.add(MassModel(), name='mass_model')

        # create_input is like a simulator input (input to the ENTIRE model)

if __name__ ==  '__main__':
    t_start             = time.time()
    # rpm                 = 5000
    # rpm_list            = np.arange(1000, 6000 + 1, 500)
    # current_list        = np.arange(50, 300 + 1, 50)
    current_list        = np.array([200])
    rpm_list            = np.array([1000])
    p                   = 12
    t                   = 0
    motor_length        = 70.e-3

    efficiency_map      = np.zeros((len(current_list), len(rpm_list)))
    input_power_map     = np.zeros((len(current_list), len(rpm_list)))
    output_power_map    = np.zeros((len(current_list), len(rpm_list)))
    torque_array        = np.zeros((len(current_list),))

    f = open('edge_deformation_data/init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    for n, current_amplitude in enumerate(current_list):
        iq                  = current_amplitude
        i_dq                = [0., iq]
        # i_abc               = [
        #     iq * np.sin(omega * t),
        #     iq * np.sin(omega * t + 2*np.pi/3),
        #     iq * np.sin(omega * t - 2*np.pi/3),
        # ]
        i_abc               = [
            iq * np.sin(0.0),
            iq * np.sin(0.0 + 2*np.pi/3),
            iq * np.sin(0.0 - 2*np.pi/3),
        ]
        Rr                  = 80.e-3
        Rs                  = 81.e-3

        fea = MotorFEA(mesh_file="mesh_files/motor_mesh_new_1", i_abc=i_abc, 
                            old_edge_coords=old_edge_coords)
        air_gap_indices         = np.arange(1,144,4)
        fea.A_z_air_gap_indices = air_gap_indices
        fea.solveMagnetostatic()

        # A_z_air_gap         = fea.A_z_air_gap
        A_z_air_gap         = fea.winding_A_z

        A_z_air_gap_shifted = np.insert(A_z_air_gap, 0, A_z_air_gap[-3:]).flatten()
        A_z_air_gap_shifted = np.delete(A_z_air_gap_shifted, [-1, -2, -3]).flatten()
        
        winding_area        = fea.winding_area_FLOAT
        magnet_area         = fea.magnet_area_FLOAT
        steel_area          = fea.steel_area_FLOAT
        print(winding_area)

        A_z_air_gap_delta   = np.array(abs(A_z_air_gap - A_z_air_gap_shifted))
        
        for m, rpm in enumerate(rpm_list):
            omega               = rpm * 2 * np.pi / 60.
            freq                = rpm * p / 120.

            # aaa = PostProcessingModel(
            #     winding_delta_A_z=A_z_air_gap_delta,
            #     motor_length=motor_length,
            #     omega=omega,
            #     winding_area=winding_area,
            #     magnet_area=magnet_area,
            #     steel_area=steel_area,
            #     phase_current_dq=i_dq,
            # )

            aaa = PostProcessingModel(
                motor_length=motor_length,
                omega=omega,
                phase_current_dq=i_dq,
                fea=fea,
                frequency=f
            )

            sim = Simulator(aaa)
            sim.run()
            print('---')
            print('wire resistance:', sim['wire_resistance'])
            print('dq voltage:',sim['phase_voltage_dq'])
            print('abc voltage:',sim['phase_voltage_abc'])
            print('dq current:',sim['phase_current_dq'])
            print('abc current:',sim['phase_current_abc'])
            print('abc flux linkage:', sim['flux_linkage_abc'])
            print('dq flux linkage:', sim['flux_linkage_dq'])
            print('input power:', sim['input_power'])
            print('output power:', sim['output_power'])
            print('output torque:', sim['output_torque'])
            print('copper_loss:', sim['copper_loss'])
            print('efficiency:', sim['efficiency'])
            print('total mass:', sim['total_mass'])
            print('---')

            efficiency_map[n, m] = sim['efficiency']
            input_power_map[n, m] = sim['input_power']
            output_power_map[n, m] = sim['output_power']

        torque_array[n] = sim['output_torque']

    # sim.visualize_implementation()
    xx, yy = np.meshgrid(rpm_list, current_list)
    contour_levels =  np.arange(70, 110  + 1, 10)

    import matplotlib.pyplot as plt
    # plt.figure(1)
    # plt.plot(motor_length_array, efficiency_array)
    # plt.ylim([0.95 * min(efficiency_array), 1.05 * max(efficiency_array)])
    # plt.xlabel('Motor Length')
    # plt.ylabel('Efficiency')

    t_end   = time.time()

    print('Runtime: ', t_end - t_start)

    plt.figure(1)
    contours = plt.contour(rpm_list, current_list, efficiency_map, levels=contour_levels)
    plt.clabel(contours,inline = True, fontsize = 8)
    plt.xlabel('RPM')
    plt.ylabel('Current Amplitude (A)')
    plt.title('Efficiency Map')

    plt.figure(2)
    contoursf = plt.contourf(rpm_list, current_list, efficiency_map, levels=np.arange(0, 120 + 1, 5))
    plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
    plt.xlabel('RPM')
    plt.ylabel('Current Amplitude (A)')
    plt.title('Efficiency Map')

    plt.figure(3)
    plt.plot(current_list, torque_array)
    plt.xlabel('Current Amplitude (A)')
    plt.ylabel('Torque (N * m)')
    plt.title('Torque')

    plt.figure(4)
    contoursf = plt.contourf(rpm_list, current_list, input_power_map, levels=np.arange(0, 8000 + 1, 500))
    plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
    plt.xlabel('RPM')
    plt.ylabel('Current Amplitude (A)')
    plt.title('Input Power Map')

    plt.figure(5)
    contoursf = plt.contourf(rpm_list, current_list, output_power_map, levels=np.arange(0, 8000 + 1, 250))
    plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
    plt.xlabel('RPM')
    plt.ylabel('Current Amplitude (A)')
    plt.title('Output Power Map')

    plt.show()
        