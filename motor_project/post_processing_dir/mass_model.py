import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

# total_mass = mass_copper + mass_steel + mass_magnet
# mass_i = density_i * area_i * length

class MassModel(Model):
    def define(self):
        p = 12
        s = p * 3
        motor_length    = self.declare_variable(name='motor_length')

        copper_density      = self.declare_variable(name='copper_density', val=8960.) # kg/m^3
        steel_density       = self.declare_variable(name='steel_density', val=7650.) # kg/m^3
        neodymium_density   = self.declare_variable(name='neodymium_density', val=7500.) # kg/m^3
        
        winding_area    = self.declare_variable(name='winding_area') # single instance of winding area
        magnet_area     = self.declare_variable(name='magnet_area') # single instance of magnet area
        steel_area      = self.declare_variable(name='steel_area') # calculated in FEniCS

        copper_area     = 2 * s * winding_area
        copper_area     = self.register_output(
            name='copper_area',
            var=copper_area
        ) # 2 * s * winding_area

        neodymium_area  = p * magnet_area
        neodymium_area  = self.register_output(
            name='neodymium_area',
            var=neodymium_area
        ) # p * magnet_area

        total_mass = motor_length * (
            copper_area * copper_density + \
            steel_area * steel_density + \
            neodymium_area * neodymium_density
        )

        total_mass = self.register_output(
            name='total_mass',
            var=total_mass
        )
        
if __name__ == '__main__':
    asdf = MassModel()
    sim = Simulator(asdf)

    sim.run()



