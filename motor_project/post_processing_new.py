import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
import time

from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.electrical_model import ElectricalModel
from post_processing_dir.power_loss_model import PowerLossModel
from post_processing_dir.mass_model import MassModel
from post_processing_dir.vector_potential_model_new import VectorPotentialModel
from post_processing_dir.torque_loss_model import TorqueLossModel
from post_processing_dir.time_average_model import TimeAverageModel

from A_z_air_gap_model import AzAirGapModel
from area_model import AreaModel
from flux_influence_ec_loss_model import FluxInfluenceECModel
from flux_influence_hysteresis_loss_model import FluxInfluenceHModel

from magnetic_states_model import MagneticStatesModel

from motor_fea import *

class PostProcessingModelDuplicate(Model):
    def define(self):
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        # final form should look like this
        # variables like winding_area will be defined as CSDL variables in the FEA
        # file, so they don't need to be declared here

class PostProcessingModel(Model):
    def initialize(self):
        self.parameters.declare('motor_length')
        self.parameters.declare('omega')
        self.parameters.declare('current_amplitude')
        self.parameters.declare('fea_list')
        self.parameters.declare('frequency')
        self.parameters.declare('angles')

    # def my_function(self, A_z, uhat, fea):
    #     e = FluxInfluenceEC(fea=self.fea)
    #     B_influence_ec = csdl.custom(A_z, uhat, op=e)
    #     return B_influence_ec

    def define(self):

        motor_length        = self.parameters['motor_length']
        omega               = self.parameters['omega']
        frequency           = self.parameters['frequency']
        current_amplitude_i = self.parameters['current_amplitude']
        fea_list            = self.parameters['fea_list']
        angles              = self.parameters['angles']
        instances           = len(angles)

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
        
        current_amplitude   = self.declare_variable(
            name='current_amplitude',
            val=current_amplitude_i,
            shape=(1,)
        )

        phase_current_dq    = self.create_output(
            name='phase_current_dq',
            val=0.,
            shape=(2,),
        )

        slot_fill_factor    = self.declare_variable(
            'slot_fill_factor',
            val=0.6
        )

        phase_current_dq[1] = current_amplitude

        flux_linkage_sign = [  # formatted as [B A C]
            [1, 0, 1],
            [1, 0, 1], # B is zero here
            [1, 0, 1],
            [1, 0, 1], # A is zero here
            [1, 0, 1],
            [1, 0, 1], # C is zero here
        ]
        
        num_windings = self.declare_variable(name='num_windings', shape=(1,), val=13.)
        
        # RU'S FEniCS MODELS
        area_model          = self.add(AreaModel(fea=fea_list[0]), name='area_model')
        
        for i in range(instances):
            # self.register_output(
            #     'flux_influence_ec_model_{}'.format(i+1),
            #     self.my_function(A_z, uhat, fea_list[i])
            # ) # alternate way of setting  up variable
            flux_influence_ec_model = self.add(FluxInfluenceECModel(fea=fea_list[i]), name='flux_influence_ec_model_{}'.format(i+1), promotes=[])
            flux_influence_h_model  = self.add(FluxInfluenceHModel(fea=fea_list[i]), name='flux_influence_h_model_{}'.format(i+1), promotes=[])
            A_z_input           = self.declare_variable(name='A_z_{}_input'.format(i+1), val=fea_list[i].A_z.vector().get_local())

            A_z_air_gap_model   = self.add(AzAirGapModel(fea=fea_list[i]), name='A_z_air_gap_model_{}'.format(i+1), promotes=[])
            vector_pot_model    = self.add(VectorPotentialModel(magnet_ref_dir=flux_linkage_sign[i]), name='vector_pot_model_{}'.format(i+1), promotes=[])
            electrical_model    = self.add(ElectricalModel(theta=angles[i], flux_linkage_sign=flux_linkage_sign[i]), name='electrical_model_{}'.format(i+1), promotes=[])

        time_average_model  = self.add(TimeAverageModel(instances=instances), name='time_average_model')

        # SET CONNECTIONS
        for i in range(instances):
            self.declare_variable('B_influence_ec_{}'.format(i+1)) 
            self.declare_variable('B_influence_hysteresis_{}'.format(i+1))

            # FLUX INFLUENCE CONNECTIONS
            self.connect('flux_influence_ec_model_{}.B_influence_ec'.format(i+1),'B_influence_ec_{}'.format(i+1))
            self.connect('flux_influence_h_model_{}.B_influence_hysteresis'.format(i+1),'B_influence_hysteresis_{}'.format(i+1))

            # ELECTRICAL MODEL 

            self.connect('phase_current_dq', 'electrical_model_{}.phase_current_dq'.format(i+1))
            # self.connect('winding_area', 'electrical_model_{}.winding_area'.format(i+1))
            # self.connect('motor_length', 'electrical_model_{}.motor_length'.format(i+1))
            # self.connect('frequency', 'electrical_model_{}.frequency'.format(i+1))
            # self.connect('slot_fill_factor', 'electrical_model_{}.slot_fill_factor'.format(i+1))
            self.connect('A_z_air_gap_model_{}.A_z_air_gap'.format(i+1), 'vector_pot_model_{}.A_z_air_gap'.format(i+1))
            self.connect('vector_pot_model_{}.flux_linkage_a_i'.format(i+1), 'electrical_model_{}.flux_linkage_a_i'.format(i+1))
            self.connect('vector_pot_model_{}.flux_linkage_b_i'.format(i+1), 'electrical_model_{}.flux_linkage_b_i'.format(i+1))
            self.connect('vector_pot_model_{}.flux_linkage_c_i'.format(i+1), 'electrical_model_{}.flux_linkage_c_i'.format(i+1))
            self.connect('electrical_model_{}.electromagnetic_torque'.format(i+1),'electromagnetic_torque_{}'.format(i+1))
            self.connect('electrical_model_{}.input_power'.format(i+1),'input_power_{}'.format(i+1))

        self.declare_variable('wire_resistance')
        self.connect('electrical_model_1.wire_resistance', 'wire_resistance')
        # SIMPLE P-P MODELS
        power_loss_model    = self.add(PowerLossModel(), name='power_loss_model')
        torque_loss_model   = self.add(TorqueLossModel(), name='torque_loss_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        mass_model          = self.add(MassModel(), name='mass_model')


if __name__ ==  '__main__':
    t_start             = time.time()
    # rpm                 = 5000
    # rpm_list            = np.arange(1000, 6000 + 1, 500)
    rpm_list            = np.array([1000])
    # current_list        = np.arange(50, 300 + 1, 50)
    current_list        = np.array([20])
    # current_list        = np.array([0.001])
    
    p                   = 12
    t                   = 0
    motor_length        = 70.e-3
    angle_shift         = 2.5
    angles              = (angle_shift + np.arange(0,30,5)) * np.pi / 180
    instances           = len(angles)
    instances           = 1

    efficiency_map      = np.zeros((len(current_list), len(rpm_list)))
    input_power_map     = np.zeros((len(current_list), len(rpm_list)))
    output_power_map    = np.zeros((len(current_list), len(rpm_list)))
    torque_array        = np.zeros((len(current_list),))
    flux_linkage_abc    = np.zeros((instances, 3))

    for n, current_amplitude in enumerate(current_list):
        print('current instance:', n)
        print('-------')
        iq                  = current_amplitude
        i_dq                = [0., iq]

        Rr                  = 80.e-3
        Rs                  = 81.e-3

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
            updateR(fea.iq, iq / 0.00016231)
            fea.solveMagnetostatic(report=True)
            

            f = open('A_z_air_gap_coords_{}.txt'.format(i+1), 'r+')
            A_z_air_gap_coords = np.fromstring(f.read(), dtype=float, sep=' ')
            f.close()
            
            fea.edge_deltas = 0.1*edge_deltas # doesn't matter right now
            fea.A_z_air_gap_indices = fea.locateAzIndices(A_z_air_gap_coords)
            print('A_z air gap values: ', fea.extractAzAirGap())
            fea_list.append(fea)

            # SAVING DATA FOR TESTING
            if False:
                vtkfile_A_z = File('post_proc_viz/Magnetic_Vector_Potential_20A_{}.pvd'.format(i+1))
                vtkfile_B = File('post_proc_viz/Magnetic_Flux_Density_20A_{}.pvd'.format(i+1))
                vtkfile_A_z << fea.A_z
                vtkfile_B << fea.B

        # end loop for instances
        # exit()
        
        for m, rpm in enumerate(rpm_list):
            print('-------')
            print('rpm instance:', m)
            print('-------')
            omega               = rpm * 2 * np.pi / 60.
            freq                = rpm * p / 120.

            aaa = PostProcessingModel(
                motor_length=motor_length,
                omega=omega,
                current_amplitude=iq,
                fea_list=fea_list,
                frequency=freq,
                angles=np.array(angles[:instances]) * p / 2
            )

            sim = Simulator(aaa)
            for i in range(instances):
                sim['A_z_air_gap_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                sim['flux_influence_ec_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                sim['flux_influence_ec_model_{}.uhat'.format(i+1)] = fea_list[i].uhat.vector().get_local()
                sim['flux_influence_h_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                sim['flux_influence_h_model_{}.uhat'.format(i+1)] = fea_list[i].uhat.vector().get_local()

            sim.run()

            if m == 0:
                for i in range(instances):
                    flux_linkage_abc[i,:] = sim['electrical_model_{}.flux_linkage_abc'.format(i+1)] 

            print('-----')
            print('flux_influence_ec_model_1:', sim['flux_influence_ec_model_1.B_influence_ec']) #  comes out as zero for some reason
            print('-----')

            print('-----')
            print('avg_flux_influence_h:', sim['avg_flux_influence_h'])
            print('flux_influence_h_list:', sim['flux_influence_h_list'])
            print('avg_flux_influence_ec:', sim['avg_flux_influence_ec'])
            print('flux_influence_ec_list:', sim['flux_influence_ec_list'])
            print('-----')
            for i in range(instances):
                print('EM torque {}'.format(i+1), sim['electromagnetic_torque_{}'.format(i+1)])
                print('Input power {}'.format(i+1), sim['electrical_model_{}.input_power'.format(i+1)])
                print('DQ flux linkage {}'.format(i+1), sim['electrical_model_{}.flux_linkage_dq'.format(i+1)])
                print('ABC flux linkage {}'.format(i+1), sim['electrical_model_{}.flux_linkage_abc'.format(i+1)])
                print('A-phase delta A_z {}'.format(i+1), sim['vector_pot_model_{}.A_z_delta_A'.format(i+1)])
                print('B-phase delta A_z {}'.format(i+1), sim['vector_pot_model_{}.A_z_delta_B'.format(i+1)])
                print('C-phase delta A_z {}'.format(i+1), sim['vector_pot_model_{}.A_z_delta_C'.format(i+1)])
                # print('ABC Current {}'.format(i+1), sim['electrical_model_{}.phase_current_abc'.format(i+1)])
                print('-----')
            print('avg EM torque:', sim['avg_electromagnetic_torque'])
            print('-----')
            print('dq current:',sim['phase_current_dq'])
            print('skin depth:', sim['wire_skin_depth'])
            print('wire radius:', sim['wire_radius'])
            print('AC resistance:', sim['wire_resistance_AC'])
            print('DC resistance:', sim['electrical_model_1.wire_resistance'])
            print('output power:', sim['output_power'])
            print('output torque:', sim['output_torque'])
            print('windage_loss:', sim['windage_loss'])
            print('hysteresis_loss:', sim['hysteresis_loss'])
            print('eddy_current_loss:', sim['eddy_current_loss'])
            print('copper_loss:', sim['copper_loss'])
            print('efficiency:', sim['efficiency'])
            print('total mass:', sim['total_mass'])
            print('---')
            print('electrical model frequency:', sim['electrical_model_1.frequency'])

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

    if False:
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

    if True:
        # for n, current_amplitude in enumerate(current_list):
        #     plt.figure(6)
        #     plt.plot(rpm_list, efficiency_map[n,:], label='Current amplitude = {} A'.format(current_amplitude))
        #     plt.xlabel('RPM')
        #     plt.ylabel('Efficiency (%)')
        #     plt.title('RPM vs. Efficiency on Current Amplitude Isolines')
        #     plt.legend()
        print(flux_linkage_abc)
        plt.figure(100)
        plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,0], 'r-*', label='Phase a')
        plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,1], 'g-*',label='Phase b')
        plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,2], 'b-*', label='Phase c')
        plt.ylabel('Flux linkage')
        plt.xlabel('Rotor mechanical angle')
        plt.title('Flux linkage of each phase (280A)')
        plt.grid()
        plt.legend()


    plt.show()
        
# ------------------------------------ NOTE: ------------------------------------
#   - in the lsdo_mesh method to output node indices, also add the coordinates as an output
