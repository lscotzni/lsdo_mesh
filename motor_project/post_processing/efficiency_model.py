import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class EfficiencyModel(Model):
    def define(self):
        omega = self.declare_variable('omega')
        winding_resistance = self.declare_variable('winding_resistance',val=1.5)
        input_power = self.declare_variable('input_power')
        output_torque =  self.declare_variable('output_torque')
        current_amplitude = self.declare_variable('current_amplitude')

        # output_power = self.create_output(
        #     name='output_power',
        #     val=output_torque * omega
        # )

        copper_loss = 3 * current_amplitude**2 * winding_resistance
        copper_loss = self.register_output(
            name='copper_loss',
            var=copper_loss
        )

        efficiency = input_power / (input_power + copper_loss)
        efficiency = self.register_output(
            name='efficiency',
            var=efficiency
        )


if __name__ == '__main__':
    aaa = EfficiencyModel()
    sim = Simulator(aaa)
    sim.run()
    print(sim['copper_loss'])
    print(sim['efficiency'])

