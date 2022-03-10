import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class TorqueLossModel(Model):
    def define(self):
        avg_electromagnetic_torque  = self.declare_variable('avg_electromagnetic_torque')
        omega                       = self.declare_variable('omega')
        hysteresis_loss             = self.declare_variable('hysteresis_loss')
        windage_loss                = self.declare_variable('windage_loss')
        eddy_current_loss           = self.declare_variable('eddy_current_loss')

        torque_loss                 = (hysteresis_loss + windage_loss + eddy_current_loss) / omega

        torque_loss     = self.register_output(
            name='torque_loss',
            var=torque_loss
        )

        output_torque               = avg_electromagnetic_torque - torque_loss

        output_torque   = self.register_output(
            name='output_torque',
            var=output_torque
        )