import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class DeltaA_zModel(Model):
    def define(self):
        winding_A_z = self.declare_variable(
            name='winding_A_z',
            shape=(72,)
        )

        winding_delta_A_z = self.create_output(
            name='winding_delta_A_z',
            shape=(36,)
        )