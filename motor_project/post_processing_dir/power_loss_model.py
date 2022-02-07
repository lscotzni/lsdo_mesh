import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class PowerLossModel(Model):
    def define(self):

        # OLD COPPER LOSS:
        # copper_loss = (current_amplitude/num_windings)**2 * wire_resistance
        # # copper_loss = 2 * (current_amplitude)**2 * wire_resistance
        # copper_loss = self.register_output(
        #     name='copper_loss',
        #     var=copper_loss
        # )
        
        # COPPER LOSS
        copper_resistivity      = self.declare_variable('copper_resistivity')
        copper_permeability     = self.declare_variable('copper_permeability')
        frequency               = self.declare_variable('frequency')
        wire_radius             = self.declare_variable('wire_radius')
        wire_resistance         = self.declare_variable('wire_resistance')
        current_amplitude       = self.declare_variable('current_amplitude')
        
        wire_skin_depth = (copper_resistivity / (np.pi * frequency * copper_permeability))**(0.5)

        wire_resistance_AC = wire_resistance / ((2*wire_skin_depth/wire_radius) - (wire_skin_depth/wire_radius)**(0.5))

        copper_loss     = (current_amplitude / np.sqrt(2))**2 * wire_resistance_AC

        copper_loss     = self.register_output(
            name='copper_loss',
            var=copper_loss
        )

        # EDDY CURRENT LOSS
        lamination_thickness    = self.declare_variable('lamination_thickness', val=0.2e-3)
        B_rms_2                 = self.declare_variable('B_rms_2')
        steel_conductivity      = self.declare_variable('steel_conductivity', val=2.17e6) # grain oriented electrical steel
        steel_area              = self.declare_variable('steel_area')
        motor_length            = self.declare_variable('motor_length')

        eddy_current_loss = np.pi**2 / 6 * steel_area * motor_length * B_rms_2**2 * frequency**2 * \
            lamination_thickness**2 * steel_conductivity

        eddy_current_loss = self.register_output(
            name='eddy_current_loss',
            var=eddy_current_loss
        )

        # HYSTERESIS LOSS
        steel_hysteresis_coeff  = self.declare_variable('steel_hysteresis_coeff', val=1.91)

        hysteresis_loss = steel_hysteresis_coeff * steel_area * motor_length * B_rms_2**2 * frequency

        hysteresis_loss = self.register_output(
            name='hysteresis_loss',
            var=hysteresis_loss
        )

        # PM LOSSES
        magnet_area         = self.declare_variable('magnet_area')
        magnet_width        = self.declare_variable('magnet_width')
        magnet_resistivity  = self.declare_variable('magnet_resistivity')

        magnet_loss     = magnet_area * motor_length * magnet_width**2 * B_rms_2**2 * frequency**2 \
            / (12 * magnet_resistivity)

        magnet_loss     = self.register_output(
            name='magnet_loss',
            var=magnet_loss
        )

        # WINDAGE LOSSES
        air_density         = self.declare_variable('air_density', val=1.204)
        air_viscosity       = self.declare_variable('air_viscosity', val=1.825e-5)
        rotor_radius        = self.declare_variable('rotor_radius', val = 80.e-3)
        air_gap_depth       = self.declare_variable('air_gap_depth', val=1.e-3)
        omega               = self.declare_variable('omega')

        azimuthal_vel       = self.register_output(
            name='azimuthal_vel',
            var=rotor_radius * omega
        )

        air_gap_Re          = air_density * azimuthal_vel * air_gap_depth / air_viscosity
        air_gap_Re          = self.register_output(
            name='air_gap_Re',
            var=air_gap_Re
        )

        friction_coeff      = 0.0152 /  (air_gap_Re)**(0.24)
        friction_coeff = self.register_output(
            name='friction_coeff',
            var=friction_coeff
        )

        windage_loss        = np.pi * air_density * friction_coeff * motor_length * omega * rotor_radius**2 # UNFINISHED
        windage_loss = self.register_output(
            name='windage_loss',
            var=windage_loss
        )

        # STRAY LOSSES
        input_power         = self.declare_variable('input_power')
        stray_loss          = input_power * 0.01

        stray_loss          = self.register_output(
            name='stray_loss',
            var=stray_loss
        )








