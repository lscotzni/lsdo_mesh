import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class EfficiencyModel(Model):
    def define(self):
        omega               = self.declare_variable('omega', shape=(1,))
        wire_resistance     = self.declare_variable('wire_resistance')
        input_power         = self.declare_variable('input_power')
        output_torque       = self.declare_variable('output_torque', shape=(1,))
        current_amplitude   = self.declare_variable('current_amplitude', 282.78/3)
        num_windings        = self.declare_variable('num_windings', 39)

        output_power = self.register_output(
            name='output_power',
            var=output_torque * omega
        )

        copper_loss = (current_amplitude/num_windings)**2 * wire_resistance
        # copper_loss = (current_amplitude/13)**2 * wire_resistance
        copper_loss = self.register_output(
            name='copper_loss',
            var=copper_loss
        )

        # efficiency = input_power / (input_power + copper_loss) # input power doesn't take into account copper loss
        # efficiency = output_power / (input_power + copper_loss)
        efficiency = output_power / (input_power + copper_loss) * 100
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

