import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class VectorPotentialModel(Model):
    def initialize(self):
        self.parameters.declare('A_z_air_gap_input')
        self.parameters.declare('magnet_ref_dir')
        self.parameters.declare('instance')

    def define(self):
        s       = 36

        # A_z_air_gap_input   = self.parameters['A_z_air_gap_input']
        magnet_ref_dir      = self.parameters['magnet_ref_dir']
        instance            = self.parameters['instance']

        A_z_air_gap = self.declare_variable(
            name='A_z_air_gap_{}'.format(instance),
            shape=(s,),
            # val=A_z_air_gap_input
        )

        A_z_delta_array = self.create_output(
            name='A_z_delta_array_{}'.format(instance),
            shape=(s,)
        )

        for i in range(s - 3):
            A_z_delta_array[i] = A_z_air_gap[i+3] - A_z_air_gap[i]

        for i in range(3):
            A_z_delta_array[s - 3 + i] = A_z_air_gap[i] - A_z_air_gap[s - 3 + i]

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

        A_z_delta_B[:]      = A_z_delta_array[0::3]
        A_z_delta_A[:]      = A_z_delta_array[1::3]
        A_z_delta_C[:]      = A_z_delta_array[2::3]

        flux_linkage_list_b     = self.create_output(
            name='flux_linkage_list_b_{}'.format(instance),
            shape=(int(s/3),)
        )

        flux_linkage_list_a     = self.create_output(
            name='flux_linkage_list_a_{}'.format(instance),
            shape=(int(s/3),)
        )

        flux_linkage_list_c     = self.create_output(
            name='flux_linkage_list_c_{}'.format(instance),
            shape=(int(s/3),)
        )

        for i in range(int(s/3)):
            flux_linkage_list_b[i] = (A_z_delta_B[i] * (-1)**(magnet_ref_dir[0] + i - 1))
            flux_linkage_list_a[i] = (A_z_delta_A[i] * (-1)**(magnet_ref_dir[1] + i - 1))
            flux_linkage_list_c[i] = (A_z_delta_C[i] * (-1)**(magnet_ref_dir[2] + i - 1))

        flux_linkage_a_i    = self.register_output(
            name='flux_linkage_a_i_{}'.format(instance),
            var=csdl.sum(flux_linkage_list_a)
        )

        flux_linkage_b_i    = self.register_output(
            name='flux_linkage_b_i_{}'.format(instance),
            var=csdl.sum(flux_linkage_list_b)
        )

        flux_linkage_c_i    = self.register_output(
            name='flux_linkage_c_i_{}'.format(instance),
            var=csdl.sum(flux_linkage_list_c)
        )

if __name__ == '__main__':
    aaa = VectorPotentialModel(
        A_z_air_gap_input=[0.1, 0.2, 0.3, -0.1, -0.2, -0.3] * 6,
        magnet_ref_dir=[0, 1, 0] 
    )
    sim = Simulator(aaa)

    sim.run()

    print(sim['A_z_air_gap'])
    print(sim['A_z_delta_array'])
    print(sim['A_z_delta_A'], sim['A_z_delta_A'].shape)
    print(sim['flux_linkage_list_a'], sim['flux_linkage_list_a'].shape)
    print(sim['flux_linkage_a_i'], sim['flux_linkage_a_i'].shape)

    print(sim['flux_linkage_list_b'], sim['flux_linkage_list_b'].shape)
    print(sim['flux_linkage_b_i'], sim['flux_linkage_b_i'].shape)

        

