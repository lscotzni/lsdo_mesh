import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class VectorPotentialModel(Model):
    def initialize(self):
        self.parameters.declare('instance')

    def define(self):
        s       = 36
        instance = self.parameters['instance']

        A_z_air_gap = self.declare_variable(
            name='A_z_air_gap_{}'.format(instance),
            shape=(s,),
            # val=np.arange(1,37)
        )

        A_z_delta_array = self.create_output(
            name='A_z_delta_array_{}'.format(instance),
            shape=(s,)
        )

        for i in range(s - 3):
            A_z_delta_array[i] = A_z_air_gap[i+3] - A_z_air_gap[i]

        for i in range(3):
            A_z_delta_array[s - 3 + i] = A_z_air_gap[i] - A_z_air_gap[s - 3 + i]

        A_z_delta_array_abs = self.create_output(
            name='A_z_delta_array_abs_{}'.format(instance),
            shape=(s,)
        )

        for i in range(s):
            A_z_delta_array_abs[i] = (((A_z_delta_array[i])**2)**(1/2))

        # A_z_delta_array_abs[:]     = A_z_delta_array[:]
            
        A_z_delta_A = self.create_output(
            name='A_z_delta_A_{}'.format(instance),
            shape=(int(s/3),)
        )

        A_z_delta_B = self.create_output(
            name='A_z_delta_B_{}'.format(instance),
            shape=(int(s/3),)
        )

        A_z_delta_C = self.create_output(
            name='A_z_delta_C_{}'.format(instance),
            shape=(int(s/3),)
        )

        A_z_delta_B[:]      = A_z_delta_array_abs[0::3]
        A_z_delta_A[:]      = A_z_delta_array_abs[1::3]
        A_z_delta_C[:]      = A_z_delta_array_abs[2::3]

        flux_linkage_a_i    = self.register_output(
            name='flux_linkage_a_i_{}'.format(instance),
            var=csdl.sum(A_z_delta_A)
        )

        flux_linkage_b_i    = self.register_output(
            name='flux_linkage_b_i_{}'.format(instance),
            var=csdl.sum(A_z_delta_B)
        )

        flux_linkage_c_i    = self.register_output(
            name='flux_linkage_c_i_{}'.format(instance),
            var=csdl.sum(A_z_delta_C)
        )

if __name__ == '__main__':
    aaa = VectorPotentialModel()
    sim = Simulator(aaa)

    sim.run()

    print(sim['A_z_delta_A'], sim['A_z_delta_A'].shape)
    print(sim['A_z_delta_A'])
    print(sim['A_z_delta_A'])
    print(sim['flux_linkage_a_i'])

        

