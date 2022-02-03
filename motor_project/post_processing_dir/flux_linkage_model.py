import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class FluxLinkageModel(Model):
    def initialize(self):
        self.parameters.declare('winding_delta_A_z')

    def define(self):
        s = 36
        # motor_length = self.declare_variable(name='motor_length', shape=(1,))
        # num_windings = self.declare_variable(name='num_windings', shape=(1,))

        # winding_delta_A_z = self.parameters['winding_delta_A_z']
        # winding_delta_A_z = self.declare_variable(
        #     name='winding_delta_A_z',
        #     val=winding_delta_A_z,
        #     shape=(s,)
        # )
        winding_delta_A_z = self.declare_variable('winding_delta_A_z', shape=(s,))

        delta_A_z_a =  self.create_output(
            name='delta_A_z_a',
            shape=(int(s/3),)
        )

        delta_A_z_b =  self.create_output(
            name='delta_A_z_b',
            shape=(int(s/3),)
        )

        delta_A_z_c =  self.create_output(
            name='delta_A_z_c',
            shape=(int(s/3),)
        )

        # for i in range(int(s/9)):
        #     delta_A_z_a[3*i:3*i+3] = winding_delta_A_z[9*i:9*i + 3]
        #     delta_A_z_b[3*i:3*i+3] = winding_delta_A_z[9*i + 3:9*i + 6]
        #     delta_A_z_c[3*i:3*i+3] = winding_delta_A_z[9*i + 6:9*i + 9]

        delta_A_z_b[:]  =  winding_delta_A_z[0::3]
        delta_A_z_a[:]  =  winding_delta_A_z[1::3]
        delta_A_z_c[:]  =  winding_delta_A_z[2::3]

        # delta_A_z_a_sum =  self.register_output(
        #     name='delta_A_z_a_sum',
        #     var=csdl.sum(delta_A_z_a),
        #     # shape=(1,)
        # )
        # delta_A_z_b_sum =  self.register_output(
        #     name='delta_A_z_b_sum',
        #     var=csdl.sum(delta_A_z_b),
        #     # shape=(1,)
        # )
        # delta_A_z_c_sum =  self.register_output(
        #     name='delta_A_z_c_sum',
        #     var=csdl.sum(delta_A_z_c),
        #     # shape=(1,)
        # )        

        # flux_linkage_a =  self.register_output(
        #     name='flux_linkage_a',
        #     var=delta_A_z_a_sum[0] * motor_length * num_windings
        # )
        # flux_linkage_b =  self.register_output(
        #     name='flux_linkage_b',
        #     var=delta_A_z_b_sum[0] * motor_length * num_windings
        # )
        # flux_linkage_c =  self.register_output(
        #     name='flux_linkage_c',
        #     var=delta_A_z_c_sum[0] * motor_length * num_windings
        # )        

        flux_linkage_a_i =  self.register_output(
            name='flux_linkage_a_i',
            var=csdl.sum(delta_A_z_a)
        )
        
        flux_linkage_b_i =  self.register_output(
            name='flux_linkage_b_i',
            var=csdl.sum(delta_A_z_b)
        )
        flux_linkage_c_i =  self.register_output(
            name='flux_linkage_c_i',
            var=csdl.sum(delta_A_z_c)
        )

        # flux_linkage_abc = self.create_output(
        #     name='flux_linkage_abc',
        #     shape=(3,)
        # )
        # print(flux_linkage_a_i.shape)
        # flux_linkage_abc[0] = flux_linkage_a_i[0] # LHS = (1,), RHS = (12,)
                
if __name__ == '__main__':
    aaa = FluxLinkageModel(winding_delta_A_z = np.arange(1,37))
    sim = Simulator(aaa)

    sim.run()

    # print(sim['winding_delta_A_z'])
    print(sim['delta_A_z_a'], sim['delta_A_z_a'].shape)
    print(sim['delta_A_z_b'])
    print(sim['delta_A_z_c'])
    print(sim['flux_linkage_a_i'].shape)

    
