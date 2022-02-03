import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class VectorPotentialModel(Model):
    def define(self):
        s       = self.declare_variable('stator_teeth')

        A_z_array = self.declare_variable(
            name='A_z_array',
            shape=(s,)
        )

        A_z_delta_array = self.create_output(
            name='A_z_delta_array',
            shape=(s,)
        )

        for i in range(s - 3):
            A_z_delta_array[i] = A_z_array[i+3] - A_z_array[i]

        for i in range(3):
            A_z_delta_array[s - 3 + i] = A_z_array[i] - A_z_array[s - 3 + i]

        A_z_delta_A = self.create_output(
            name='A_z_delta_A',
            shape=(s/3,)
        )

        A_z_delta_B = self.create_output(
            name='A_z_delta_A',
            shape=(s/3,)
        )

        A_z_delta_C = self.create_output(
            name='A_z_delta_A',
            shape=(s/3,)
        )

        A_z_delta_B[:]      = A_z_delta_array[0::3]
        A_z_delta_A[:]      = A_z_delta_array[1::3]
        A_z_delta_C[:]      = A_z_delta_array[2::3]
        
        # ===============================================================
        # NOTE: need to implement absolute values for the summations
        # can do something like this:
        # for i in range(s/3):
        #   if A_z_delta_B[i] < 0:
        #       A_z_delta_B[i] = -1 * A_z_delta_B[i]
        #   repeat for all phases
        # ===============================================================

        flux_linkage_a_i    = self.register_output(
            name='flux_linkage_a_i',
            var=csdl.sum(A_z_delta_A)
        )

        flux_linkage_b_i    = self.register_output(
            name='flux_linkage_b_i',
            var=csdl.sum(A_z_delta_B)
        )

        flux_linkage_c_i    = self.register_output(
            name='flux_linkage_c_i',
            var=csdl.sum(A_z_delta_C)
        )


        

