import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

# total_mass = mass_copper + mass_steel + mass_magnet
# mass_i = density_i * area_i * length

class MassModel(Model):
    def define(self):
        motor_length    = self.declare_variable(name='motor_length')

        copper_density  = self.declare_variable(name='copper_density', val=8960.) # kg/m^3
        steel_density   = self.declare_variable(name='steel_density', val=7650.) # kg/m^3
        magnet_density  = self.declare_variable(name='magnet_density', val=7500.) # kg/m^3
        
        winding_area    = self.declare_variable(name='winding_area') # single instance of winding area
        magnet_area     = self.declare_variable(name='magnet_area') # single instance of magnet area
        steel_area      = self.declare_variable(name='steel_area') # calculated in FEniCS

        copper_area     = self.create_output(name='copper_area') # 2 * s * winding_area
        neodymium_area  = self.create_output(name='neodymium_area') # p * magnet_area
        
        


