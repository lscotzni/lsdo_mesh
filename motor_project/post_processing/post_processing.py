import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

from efficiency_model import EfficiencyModel
from flux_linkage_model import FluxLinkageModel
from electrical_model import ElectricalModel


class PostProcessingModel(Model):
    def initialize(self):
        self.parameters.declare('omega')
    def define(self):
        # omega = self.parameters['omega']
        # omega = self.declare_variable('omega', 1000)
        # omega =  self.create_output('omega', val=omega)
        omega = self.declare_variable('omega')
        
        flux_linkage_model  = self.add(FluxLinkageModel(), name='flux_linkage_model')
        electrical_model    = self.add(ElectricalModel(), name='electrical_model')
        efficiency_model    = self.add(EfficiencyModel(), name='efficiency_model')

        self.connect('omega','electrical_model.omega')
        self.connect('omega','efficiency_model.omega')
        # self.connect('electrical_model.omega','efficiency_model.omega')
        self.connect('efficiency_model.input_power', 'electrical_model.input_power')
        self.connect('efficiency_model.output_torque', 'electrical_model.output_torque')

if __name__ ==  '__main__':
    aaa = PostProcessingModel(omega='omega')
    sim = Simulator(aaa)
    sim.run()
        