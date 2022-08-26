import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

from motor_mesh import MotorMeshGenerator
from post_processing_dir.efficiency_model import EfficiencyModel
from post_processing_dir.flux_linkage_model import FluxLinkageModel
from post_processing_dir.electrical_model import ElectricalModel
from post_processing_dir.mass_model import MassModel

class MotorModel(Model):
    def initialize(self):
        self.parameters.declare('mesh_model_instances')

    def define(self):

        mesh_model_instances    = self.parameters['mesh_model_instances'] # holds mesh models for each rotation instance

        # MAGNET DV
        magnet_thickness        = self.create_input('magnet_thickness', shape=(1,))
        self.add_design_variable('magnet_thickness')

        magnet_angular_sweep    = self.create_input('magnet_angular_sweep', shape=(1,))
        self.add_design_variable('magnet_angular_sweep')

        # ROTOR DV
        outer_rotor_radius      = self.create_input('outer_rotor_radius', shape=(1,))
        self.add_design_variable('outer_rotor_radius')

        # STATOR DV
        outer_stator_radius     = self.create_input('outer_stator_radius', shape=(1,))
        self.add_design_variable('outer_stator_radius')

        inner_stator_radius     = self.create_input('inner_stator_radius', shape=(1,))
        self.add_design_variable('inner_stator_radius')

        tooth_angular_sweep     = self.create_input('tooth_angular_sweep', shape=(1,))
        self.add_design_variable('tooth_angular_sweep')

        tooth_height            = self.create_input('tooth_height', shape=(1,))
        self.add_design_variable('tooth_height')

        # WINDING DV
        winding_angular_sweep   = self.create_input('winding_angular_sweep', shape=(1,))
        self.add_design_variable('winding_angular_sweep')





        


        for i in range(len(mesh_model_instances)):
            self.add(mesh_model_instances[i], name='mesh_model_{}'.format(i+1), promotes=[])
            self.add(EfficiencyModel(), name='efficiency_model_{}'.format(i+1), promotes=[])
            self.add(ElectricalModel(), name='electrical_model{}'.format(i+1), promotes=[])
            self.add(FluxLinkageModel(), name='flux_linkage_model{}'.format(i+1), promotes=[])



