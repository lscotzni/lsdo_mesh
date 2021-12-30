import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.flux_linkage_model import FluxLinkageModel
from post_processing_dir.electrical_model import ElectricalModel
from post_processing_dir.mass_model import MassModel

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

    def define(self):
        winding_delta_A_z   = self.parameters['winding_delta_A_z']
        winding_area        = self.parameters['winding_area']
        motor_length        = self.parameters['motor_length']
        omega               = self.parameters['omega']
        magnet_area         = self.parameters['magnet_area']
        steel_area          = self.parameters['steel_area']
        phase_current_dq    = self.parameters['phase_current_dq']

        winding_delta_A_z   = self.declare_variable(
            name='winding_delta_A_z',
            val=winding_delta_A_z,
            shape=(36,)
        )
        winding_area        = self.declare_variable(
            name='winding_area',
            val=winding_area
        )
        motor_length        = self.declare_variable(
            name='motor_length',
            val=motor_length
        )
        omega               = self.declare_variable(
            name='omega',
            val=omega
        )
        magnet_area         = self.declare_variable(
            name='magnet_area',
            val=magnet_area
        )
        steel_area          = self.declare_variable(
            name='steel_area',
            val=steel_area
        )
        phase_current_dq    = self.declare_variable(
            name='phase_current_dq',
            val=phase_current_dq,
            shape=(2,),
        )
        
        flux_linkage_model  = self.add(FluxLinkageModel(), name='flux_linkage_model')
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        mass_model          = self.add(MassModel(), name='mass_model')

        # create_input is like a simulator input (input to the ENTIRE model)

if __name__ ==  '__main__':
    rpm                 = 1000
    omega               = rpm * 2 * np.pi / 60.
    t                   = 0

    current_amplitude   = 282.8
    iq                  = -current_amplitude / 3
    i_dq                = [0., iq]
    i_abc               = [
        -iq * np.sin(omega * t + np.pi/2),
        -iq * np.sin(omega * t - 2*np.pi/3 + np.pi/2),
        -iq * np.sin(omega * t + 2*np.pi/3 + np.pi/2),
    ]

    fea = MagnetostaticProblem(i_abc=i_abc)
    fea.solveMagnetostatic()
    winding_delta_A_z   = fea.winding_delta_A_z
    winding_area        = fea.winding_area
    magnet_area         = fea.magnet_area
    steel_area          = fea.steel_area

    # motor_length_array = np.linspace(2,10,20)
    motor_length_array = np.linspace(0.5,1.5,10)
    efficiency_array, mass_array = [], []
    motor_length = motor_length_array[0]
    for length in motor_length_array:
        aaa = PostProcessingModel(
            winding_delta_A_z=winding_delta_A_z,
            motor_length=length,
            omega=omega,
            winding_area=winding_area,
            magnet_area=magnet_area,
            steel_area=steel_area,
            phase_current_dq=i_dq,
        )
              
        sim = Simulator(aaa)
        sim.run()
        print('---')
        # print(sim['winding_delta_A_z'])
        # print(sim['flux_linkage_a_i'])
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
        efficiency_array.extend(sim['efficiency'])
        mass_array.extend(sim['total_mass'])

        # print(sim['transform_matrix_forward'])
        # print(sim['transform_matrix_reverse'])
        # print(sim['winding_delta_A_z'])
        # exit()
    
    print(efficiency_array)
    print(mass_array)


    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(motor_length_array, efficiency_array)
    plt.ylim([0.95 * min(efficiency_array), 1.05 * max(efficiency_array)])
    plt.xlabel('Motor Length')
    plt.ylabel('Efficiency')

    plt.figure(2)
    plt.plot(motor_length_array, mass_array)
    plt.ylim([0.95 * min(mass_array), 1.05 * max(mass_array)])
    plt.xlabel('Motor Length')
    plt.ylabel('Motor Mass (kg)')

    plt.show()
        