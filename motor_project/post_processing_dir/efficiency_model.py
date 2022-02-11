import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class EfficiencyModel(Model):
    def define(self):
        omega               = self.declare_variable('omega', shape=(1,))
        input_power         = self.declare_variable('input_power')
        output_torque       = self.declare_variable('output_torque', shape=(1,))
        current_amplitude   = self.declare_variable('current_amplitude')
        num_windings        = self.declare_variable('num_windings')

        # INDIVIDUAL POWER LOSSES
        copper_loss         = self.declare_variable(name='copper_loss')
        eddy_current_loss   = self.declare_variable(name='eddy_current_loss')
        hysteresis_loss     = self.declare_variable(name='hysteresis_loss')
        magnet_loss         = self.declare_variable(name='magnet_loss')
        windage_loss        = self.declare_variable(name='windage_loss')
        stray_loss          = self.declare_variable(name='stray_loss')

        output_power = self.register_output(
            name='output_power',
            var=output_torque * omega
        )

        total_power_loss = self.register_output(
            name='total_power_loss',
            var=copper_loss + eddy_current_loss + hysteresis_loss + magnet_loss + windage_loss + stray_loss
        )

        # NEED TO ADD COMPUTATIONS FOR LOAD TORQUE

        # efficiency = input_power / (input_power + copper_loss) # input power doesn't take into account copper loss
        # efficiency = output_power / (input_power + copper_loss)
        efficiency = output_power / (output_power + total_power_loss) * 100
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

