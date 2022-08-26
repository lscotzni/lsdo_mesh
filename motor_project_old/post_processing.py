import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
import time

from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.electrical_model import ElectricalModel
from post_processing_dir.power_loss_model import PowerLossModel
from post_processing_dir.mass_model import MassModel
from post_processing_dir.flux_linkage_model import FluxLinkageModel
from post_processing_dir.torque_loss_model import TorqueLossModel
from post_processing_dir.time_average_model import TimeAverageModel

from csdl_fenics_models.area_model import AreaModel
from magnetic_states_model import MagneticStatesModel

from csdl_fenics_models.A_z_air_gap_model import AzAirGapModel
from csdl_fenics_models.flux_influence_ec_loss_model import FluxInfluenceECModel
from csdl_fenics_models.flux_influence_hysteresis_loss_model import FluxInfluenceHModel

# # IMPORTS BELOW ARE CustomExplicitOperation FOR MANUAL FUNCTIONS
# from A_z_air_gap_model import AzAirGap
# from flux_influence_ec_loss_model import FluxInfluenceEC
# from flux_influence_hysteresis_loss_model import FluxInfluence


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

    def flux_influence_EC_function(self, fea):
        # A_z     = fea.A_z.vector().get_local()
        # uhat    = fea.uhat.vector().get_local()
        e  = FluxInfluenceEC(fea=fea)
        B_influence_ec = csdl.custom(A_z, uhat, op=e)
        return B_influence_ec

    def flux_influence_H_function(self, A_z, uhat, fea):
        e  = FluxInfluenceH(fea=fea)
        B_influence_h = csdl.custom(A_z, uhat, op=e)
        return B_influence_h

    def air_gap_A_z_function(self, A_z, fea):
        e  = AzAirGap(fea=self.fea)
        A_z_air_gap = csdl.custom(A_z, op=e)
        return A_z_air_gap
    
    # def my_function(self, A_z, uhat, fea):
    #     e = FluxInfluenceEC(fea=self.fea)
    #     B_influence_ec = csdl.custom(A_z, uhat, op=e)
    #     return B_influence_ec

    def define(self):

        motor_length_i      = self.parameters['motor_length']
        omega               = self.parameters['omega']
        frequency_i         = self.parameters['frequency']
        current_amplitude_i = self.parameters['current_amplitude']
        fea_list            = self.parameters['fea_list']
        angles              = self.parameters['angles']
        instances           = len(angles)

        motor_length        = self.declare_variable('motor_length', motor_length_i)
        omega               = self.declare_variable('omega', omega)
        frequency           = self.declare_variable('frequency', frequency_i)
        current_amplitude   = self.declare_variable('current_amplitude', current_amplitude_i)
        phase_current_dq    = self.create_output('phase_current_dq', val=0., shape=(2,))
        phase_current_dq[1] = current_amplitude

        # NOT NEEDED ANYMORE
        flux_linkage_sign = [  # formatted as [B A C]
            [1, 0, 1],
            [1, 0, 1], # B is zero here
            [1, 0, 1],
            [1, 0, 1], # A is zero here
            [1, 0, 1],
            [1, 0, 1], # C is zero here
        ]
        
        for i in range(instances):
            # self.register_output(
            #     'flux_influence_ec_{}'.format(i+1),
            #     self.flux_influence_EC_function(A_z, uhat, fea_list[i])
            # ) # alternate way of setting  up variable

            instance_model = Model()
            
            instance_model.declare_variable('motor_length', motor_length_i)
            instance_model.declare_variable('frequency', frequency_i)
            instance_model.declare_variable('num_windings', 13.)
            instance_model.declare_variable('slot_fill_factor', val=0.6) # only kept here bc this can potentially change
            

            instance_model.add(FluxInfluenceECModel(fea=fea_list[i]), name='flux_influence_ec_model_{}'.format(i+1))
            instance_model.add(FluxInfluenceHModel(fea=fea_list[i]), name='flux_influence_h_model_{}'.format(i+1))
            instance_model.add(AzAirGapModel(fea=fea_list[i]), name='A_z_air_gap_model_{}'.format(i+1))
            instance_model.add(AreaModel(fea=fea_list[i]), name='area_model_{}'.format(i+1))
            instance_model.add(FluxLinkageModel(), name='flux_linkage_model_{}'.format(i+1))
            instance_model.add(ElectricalModel(theta=angles[i]), name='electrical_model_{}'.format(i+1))

            self.add(instance_model, 'instance_model_{}'.format(i+1), promotes=[])

        time_average_model  = self.add(TimeAverageModel(instances=instances), name='time_average_model')
        wire_resistance     = self.declare_variable('wire_resistance')
        

        # SET CONNECTIONS
        for i in range(instances):
            self.connect('phase_current_dq', 'instance_model_{}.phase_current_dq'.format(i+1))

            # FLUX INFLUENCE CONNECTIONS
            self.connect('instance_model_{}.B_influence_ec'.format(i+1), 'B_influence_ec_{}'.format(i+1))
            self.connect('instance_model_{}.B_influence_h'.format(i+1), 'B_influence_h_{}'.format(i+1))
            self.connect('instance_model_{}.electromagnetic_torque'.format(i+1), 'electromagnetic_torque_{}'.format(i+1))
            self.connect('instance_model_{}.input_power'.format(i+1), 'input_power_{}'.format(i+1))

        self.connect('instance_model_1.wire_resistance', 'wire_resistance')

        # SIMPLE P-P MODELS
        power_loss_model    = self.add(PowerLossModel(), name='power_loss_model')
        torque_loss_model   = self.add(TorqueLossModel(), name='torque_loss_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        mass_model          = self.add(MassModel(), name='mass_model')


if __name__ ==  '__main__':
    slot_fill_factor    = 0.6
    num_windings        = 13
    t_start             = time.time()
    rpm_list            = np.arange(1000, 6000 + 1, 500)
    # rpm_list            = np.array([1000., 2000.])
    # current_list        = np.arange(20., 220. + 1., 20)
    current_list        = np.array([0.])
    
    p                   = 12
    t                   = 0
    motor_length        = 70.e-3
    angle_shift         = 2.5
    angles              = (angle_shift + np.arange(0,30,5)) * np.pi / 180
    # instances           = len(angles)
    instances           = 2

    efficiency_map      = np.zeros((len(current_list), len(rpm_list)))
    input_power_map     = np.zeros((len(current_list), len(rpm_list)))
    output_power_map    = np.zeros((len(current_list), len(rpm_list)))
    torque_array        = np.zeros((len(current_list),))

    eddy_current_loss_map   = np.zeros((len(current_list), len(rpm_list)))
    hysteresis_loss_map     = np.zeros((len(current_list), len(rpm_list)))
    copper_loss_array       = np.zeros((len(current_list),))

    flux_linkage_abc    = np.zeros((instances, 3))

    windage_loss_array      = np.zeros((len(rpm_list), ))
    windage_loss_vel_array  = np.zeros((len(rpm_list), ))
    air_gap_Ta_array        = np.zeros((len(rpm_list), ))
    air_gap_Re_array        = np.zeros((len(rpm_list), ))
    
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
            updateR(fea.iq, iq / (0.00016231 / num_windings))
            # fea.iq = iq / (0.00016231 / num_windings)
            fea.solveMagnetostatic(report=True)
            
            f = open('edge_deformation_data/A_z_air_gap_coords_{}.txt'.format(i+1), 'r+')
            A_z_air_gap_coords = np.fromstring(f.read(), dtype=float, sep=' ')
            f.close()
            
            fea.edge_deltas = 0.1*edge_deltas # doesn't matter right now
            fea.A_z_air_gap_indices = fea.locateAzIndices(A_z_air_gap_coords)
            print('A_z air gap values: ', fea.extractAzAirGap())
            fea_list.append(fea)

            # SAVING DATA FOR TESTING
            if True:
                vtkfile_A_z = File('post_proc_viz/Magnetic_Vector_Potential_{}_{}_0726TEST.pvd'.format(str(current_amplitude),i+1))
                vtkfile_B = File('post_proc_viz/Magnetic_Flux_Density_{}_{}_0726TEST.pvd'.format(str(current_amplitude),i+1))
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

            # aaa.visualize_sparsity()

            sim = Simulator(aaa)
            for i in range(instances):
                # sim['A_z_air_gap_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                # sim['flux_influence_ec_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                # sim['flux_influence_ec_model_{}.uhat'.format(i+1)] = fea_list[i].uhat.vector().get_local()
                # sim['flux_influence_h_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                # sim['flux_influence_h_model_{}.uhat'.format(i+1)] = fea_list[i].uhat.vector().get_local()

                sim['instance_model_{}.A_z'.format(i+1)] = fea_list[i].A_z.vector().get_local()
                sim['instance_model_{}.uhat'.format(i+1)] = fea_list[i].uhat.vector().get_local()

            sim.run()

            if m == 0:
                for i in range(instances):
                    flux_linkage_abc[i,:] = sim['instance_model_{}.flux_linkage_abc'.format(i+1)] 

            print('-----')
            print('avg_flux_influence_h:', sim['avg_flux_influence_h'])
            print('flux_influence_h_list:', sim['flux_influence_h_list'])
            print('avg_flux_influence_ec:', sim['avg_flux_influence_ec'])
            print('flux_influence_ec_list:', sim['flux_influence_ec_list'])
            print('-----')
            for i in range(instances):
                print('EM torque {}'.format(i+1), sim['instance_model_{}.electromagnetic_torque'.format(i+1)])
                print('Input power {}'.format(i+1), sim['instance_model_{}.input_power'.format(i+1)])
                print('DQ flux linkage {}'.format(i+1), sim['instance_model_{}.flux_linkage_dq'.format(i+1)])
                print('ABC flux linkage {}'.format(i+1), sim['instance_model_{}.flux_linkage_abc'.format(i+1)])
                print('A-phase delta A_z {}'.format(i+1), sim['instance_model_{}.A_z_delta_A'.format(i+1)])
                print('B-phase delta A_z {}'.format(i+1), sim['instance_model_{}.A_z_delta_B'.format(i+1)])
                print('C-phase delta A_z {}'.format(i+1), sim['instance_model_{}.A_z_delta_C'.format(i+1)])
                print('v-pot model flux linkage a: ', sim['instance_model_{}.flux_linkage_a_i'.format(i+1)])
                print('ABC voltage {}'.format(i+1), sim['instance_model_{}.phase_voltage_abc'.format(i+1)])

                print('-----')
            print('avg EM torque:', sim['avg_electromagnetic_torque'])
            print('-----')
            print('frequency:', sim['frequency'])
            print('dq current:',sim['phase_current_dq'])
            print('skin depth:', sim['wire_skin_depth'])
            print('wire radius:', sim['wire_radius'])
            print('AC resistance:', sim['wire_resistance_AC'])
            print('DC resistance:', sim['wire_resistance'])
            print('input power:', sim['output_power'])
            print('output power:', sim['output_power'])
            print('output torque:', sim['output_torque'])
            print('windage_loss:', sim['windage_loss'])
            print('Air gap Re: ', sim['air_gap_Re'])
            print('hysteresis_loss:', sim['hysteresis_loss'])
            print('eddy_current_loss:', sim['eddy_current_loss'])
            print('copper_loss:', sim['copper_loss'])
            print('total power loss: ', sim['total_power_loss'])
            print('azimuthal velocity: ', sim['azimuthal_vel'])
            print('efficiency:', sim['efficiency'])
            print('total mass:', sim['total_mass'])
            print('winding area: ', sim['winding_area'])
            print('steel area: ', sim['steel_area'])
            print('magnet area: ', sim['magnet_area'])

            print(' ---- Windage Loss Analysis ---- ')
            print('windage_loss:', sim['windage_loss'])
            print('Air gap Re: ', sim['air_gap_Re'])
            print('azimuthal velocity: ', sim['azimuthal_vel'])
            print('friction_coeff: ', sim['friction_coeff'])
            print('omega: ', sim['omega'])

            efficiency_map[n, m] = sim['efficiency']
            input_power_map[n, m] = sim['avg_input_power']
            output_power_map[n, m] = sim['output_power']
            eddy_current_loss_map[n,m] = sim['eddy_current_loss']
            hysteresis_loss_map[n,m] = sim['hysteresis_loss']
            if i == 0 and n == 0:
                windage_loss_array[m] = sim['windage_loss']
                windage_loss_vel_array[m] = sim['azimuthal_vel']
                air_gap_Ta_array[m]    = sim['air_gap_Ta']
                air_gap_Re_array[m]    = sim['air_gap_Re']

            
        copper_loss_array[n] = sim['copper_loss']
        torque_array[n] = sim['output_torque']

    xx, yy = np.meshgrid(rpm_list, current_list)
    contour_levels =  np.arange(60, 100  + 1, 5)

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

        # plt.figure(1)
        # contours = plt.contour(rpm_list, current_list, efficiency_map)
        # plt.clabel(contours,inline = True, fontsize = 8)
        # plt.xlabel('RPM')
        # plt.ylabel('Current Amplitude (A)')
        # plt.title('Efficiency Map')

        plt.figure(2)
        contoursf = plt.contourf(rpm_list, current_list, efficiency_map, levels=np.arange(0, 100 + 1, 5))
        plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
        plt.xlabel('RPM')
        plt.ylabel('Current Amplitude (A)')
        plt.title('Efficiency Map')

        plt.figure(3)
        contoursf = plt.contourf(rpm_list, torque_array , efficiency_map, levels=np.arange(85, 100 + 1, 1))
        contours = plt.contour(rpm_list, torque_array , efficiency_map, levels=np.arange(85, 100 + 1, 1))
        # plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
        plt.clabel(contours,inline = True, fontsize = 8, colors='black')
        plt.xlabel('RPM')
        plt.ylabel('EM Torque (Nm)')
        plt.title('Efficiency Map')

        # plt.figure(3)
        # plt.plot(current_list, torque_array)
        # plt.xlabel('Current Amplitude (A)')
        # plt.ylabel('Torque (N * m)')
        # plt.title('Torque')

        # plt.figure(4)
        # contoursf = plt.contourf(rpm_list, current_list, input_power_map, levels=np.arange(0, 8000 + 1, 500))
        # plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
        # plt.xlabel('RPM')
        # plt.ylabel('Current Amplitude (A)')
        # plt.title('Input Power Map')

        # plt.figure(5)
        # contoursf = plt.contourf(rpm_list, current_list, output_power_map, levels=np.arange(0, 8000 + 1, 250))
        # plt.clabel(contoursf,inline = True, fontsize = 8, colors='black')
        # plt.xlabel('RPM')
        # plt.ylabel('Current Amplitude (A)')
        # plt.title('Output Power Map')

    if False:
        for n, current_amplitude in enumerate(current_list):
            plt.figure(6)
            plt.plot(rpm_list, efficiency_map[n,:], '*-', linewidth=1, markersize=5, label='{} A'.format(current_amplitude))
            plt.xlabel('RPM')
            plt.ylabel('Efficiency (%)')
            plt.title('RPM vs. Efficiency (Current Amplitude Isolines)')
            plt.legend()
            plt.grid(visible=True)

            plt.figure(7)
            plt.plot(rpm_list, input_power_map[n,:] / 1000, '*-', linewidth=1, markersize=5, label='{} A'.format(current_amplitude))
            plt.xlabel('RPM')
            plt.ylabel('Input Power (kW)')
            plt.title('RPM vs. Input Power (Current Amplitude Isolines)')
            plt.legend()
            plt.grid(visible=True)

            plt.figure(8)
            plt.plot(rpm_list, output_power_map[n,:] / 1000, '*-', linewidth=1, markersize=5, label='{} A'.format(current_amplitude))
            plt.xlabel('RPM')
            plt.ylabel('Output Power (kW)')
            plt.title('RPM vs. Output Power (Current Amplitude Isolines)')
            plt.legend()
            plt.grid(visible=True)

            plt.figure(9)
            plt.plot(rpm_list, eddy_current_loss_map[n,:], '*-', linewidth=1, markersize=5, label='{} A'.format(current_amplitude))
            plt.xlabel('RPM')
            plt.ylabel('Eddy Current Loss (W)')
            plt.title('RPM vs. Eddy Current Loss (Current Amplitude Isolines)')
            plt.legend()
            plt.grid(visible=True)

            plt.figure(10)
            plt.plot(rpm_list, hysteresis_loss_map[n,:], '*-', linewidth=1, markersize=5, label='{} A'.format(current_amplitude))
            plt.xlabel('RPM')
            plt.ylabel('Hysteresis Loss (W)')
            plt.title('RPM vs. Hysteresis Loss (Current Amplitude Isolines)')
            plt.legend()
            plt.grid(visible=True)

            plt.figure(11)
            plt.plot(rpm_list, windage_loss_array, 'k', linewidth=3)
            plt.xlabel('RPM')
            plt.ylabel('Windage Loss (W)')
            plt.title('RPM vs. Windage Loss (Current Amplitude Isolines)')
            plt.grid(visible=True)

            plt.figure(12)
            plt.plot(current_list, copper_loss_array, 'k', linewidth=3)
            plt.xlabel('Current Amplitude (A)')
            plt.ylabel('Copper Loss (W)')
            plt.title('Current vs. Copper Loss')
            plt.grid(visible=True)

            plt.figure(13)
            plt.plot(current_list, torque_array, 'k', linewidth=3)
            plt.xlabel('Current Amplitude (A)')
            plt.ylabel('EM Torque  (Nm)')
            plt.title('Current vs. Electromagnetic Torque')
            plt.grid(visible=True)

            plt.figure(14)
            plt.plot(rpm_list * 2*np.pi/60 * 80e-3, windage_loss_vel_array, 'k', linewidth=3)
            plt.xlabel('Rotor edge linear velocity (m/s)')
            plt.ylabel('Average azimuthal velocity (m/s)')
            plt.title('Rotor Speed vs. Avg. Azimuthal velocity')
            plt.grid(visible=True)

            # plt.figure(15)
            # plt.plot(air_gap_Ta_array, windage_loss_array, 'k', linewidth=3)
            # plt.xlabel('Air gap Taylor Number')
            # plt.ylabel('Windage Loss (W)')
            # plt.title('Ta vs. Windage Loss')
            # plt.grid(visible=True)

            # plt.figure(16)
            # plt.plot(air_gap_Re_array, windage_loss_array, 'k', linewidth=3)
            # plt.xlabel('Air gap Reynolds Number')
            # plt.ylabel('Windage Loss (W)')
            # plt.title('Re vs. Windage Loss')
            # plt.grid(visible=True)

            print(flux_linkage_abc)
            plt.figure(100)
            plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,0], 'r-*', label='Phase a')
            plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,1], 'g-*',label='Phase b')
            plt.plot(angles[:instances] * 180 / np.pi, flux_linkage_abc[:,2], 'b-*', label='Phase c')
            plt.ylabel('Flux linkage')
            plt.xlabel('Rotor mechanical angle')
            plt.title('Flux linkage of each phase ({}A)'.format(current_amplitude))
            plt.grid()
            plt.legend()


    # plt.show()
        
# ------------------------------------ NOTE: ------------------------------------
#   - in the lsdo_mesh method to output node indices, also add the coordinates as an output
