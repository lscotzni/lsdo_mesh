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
    def define(self):
        # omega = self.declare_variable('omega', 1000)
        # omega =  self.create_output('omega', val=omega)
        winding_delta_A_z   = self.parameters['winding_delta_A_z']
        winding_area        = self.parameters['winding_area']
        motor_length        = self.parameters['motor_length']

        omega               = self.declare_variable('omega', val=2500 * 2 * np.pi / 60) # use create_input instead
        # omega               = self.create_input('omega', val=2500 * 2 * np.pi / 60)
        winding_area        = self.declare_variable('winding_area', winding_area)
        motor_length        = self.declare_variable('motor_length', motor_length)
        winding_delta_A_z   = self.declare_variable('winding_delta_A_z', winding_delta_A_z, shape=(36,))
        
        # flux_linkage_model  = self.add(FluxLinkageModel(winding_delta_A_z=winding_delta_A_z), name='flux_linkage_model')
        flux_linkage_model  = self.add(FluxLinkageModel(), name='flux_linkage_model')
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')
        mass_model          = self.add(MassModel(), name='mass_model')

        # self.connect('omega','electrical_model.omega')
        # self.connect('omega','efficiency_model.omega')
        # self.connect('electrical_model.omega','efficiency_model.omega')
        # self.connect('efficiency_model.input_power', 'electrical_model.input_power')
        # self.connect('efficiency_model.output_torque', 'electrical_model.output_torque')

        # create_input is like a simulator input (input to the ENTIRE model)

if __name__ ==  '__main__':
    import matplotlib.pyplot as plt

    fea = MagnetostaticProblem()
    fea.solveMagnetostatic()
    winding_delta_A_z   = fea.winding_delta_A_z
    winding_area        = fea.winding_area
    magnet_area         = fea.magnet_area
    steel_area          = fea.steel_area
    print(winding_area)
    # print(winding_delta_A_z.shape)
    motor_length_array = np.linspace(2,10,20)
    efficiency_array = []
    motor_length = motor_length_array[0]
    for length in motor_length_array:
        aaa = PostProcessingModel(
            winding_delta_A_z=winding_delta_A_z,
            winding_area=winding_area,
            motor_length=length,
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
        print('---')
        efficiency_array.extend(sim['efficiency'])

        print(sim['transform_matrix_forward'])
        print(sim['transform_matrix_reverse'])
        print(sim['winding_delta_A_z'])
        # exit()
    
    print(efficiency_array)

    plt.figure(1)
    plt.plot(motor_length_array, efficiency_array)
    plt.show()
        