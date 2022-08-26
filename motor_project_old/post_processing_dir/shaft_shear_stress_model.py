import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

'''
NOTE: THIS CLASS IS << NOT >> POST-PROCESSING
    - WE KNOW THE NECESSARY TORQUE AND MATERIAL SHEAR STRESS
    - SHAFT DIAMETER CAN BE PRESCRIBED AS A FUNCTION OF THE TWO
'''

class ShaftShearStressModel(Model):
    def define(self):
        method = 'diameter' # or 'shear_stress'

        shaft_shear_stress = self.declare_variable(
            name='shaft_shear_stress',
            val=1.e6
        ) # DEPENDENT ON MATERIAL OF SHAFT

        output_torque = self.declare_variable(
            name='output_torque'
        ) # COMES FROM ELECTRICAL MODEL

        min_shaft_diameter = 1.72 * (output_torque / shaft_shear_stress) ** 0.33
        min_shaft_diameter = self.register_output(
            name='min_shaft_diameter',
            var=min_shaft_diameter
        )
        # looking at minimum shaft diameter



if __name__ == '__main__':
    asdf = ShaftShearStressModel()
    sim = Simulator(asdf)

    sim.run()



# NOTES:
# D^3 = 16/pi * T / tau 
# tau = 16 / pi * T / D^3
# can set such that the RHS cannot exceed the max 
# shear stress of given material

