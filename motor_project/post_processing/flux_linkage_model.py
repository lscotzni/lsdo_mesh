import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class FluxLinkageModel(Model):
    def define(self):
        s = 36
        motor_length = self.create_input(name='motor_length')
        num_windings = self.create_input(name='num_windings')

        winding_delta_A_z = self.declare_variable(
            name='winding_delta_A_z',
            val=np.arange(1,37),
            shape=(s,)
        )
        # THIS HOLDS THE AVERAGE MAGNETIC VECTOR POTENTIAL VALUE IN A WINDING
        # WE NEED TO TAKE THE DELTA

        delta_A_z_a =  self.create_output(
            name='delta_A_z_a',
            shape=(int(s/3),)
        )
        # delta_A_z_a[:] = winding_delta_A_z[:12]

        delta_A_z_b =  self.create_output(
            name='delta_A_z_b',
            shape=(int(s/3),)
        )

        delta_A_z_c =  self.create_output(
            name='delta_A_z_c',
            shape=(int(s/3),)
        )

        for i in range(int(s/9)):
            delta_A_z_a[3*i:3*i+3] = winding_delta_A_z[9*i:9*i + 3]
            delta_A_z_b[3*i:3*i+3] = winding_delta_A_z[9*i + 3:9*i + 6]
            delta_A_z_c[3*i:3*i+3] = winding_delta_A_z[9*i + 6:9*i + 9]

        # flux_linkage_a_list =  self.create_output(
        #     name='flux_linkage_a_list',
        #     shape=(int(s/3),)
        # )

        # flux_linkage_b_list =  self.create_output(
        #     name='flux_linkage_b_list',
        #     shape=(int(s/3),)
        # )

        # flux_linkage_c_list =  self.create_output(
        #     name='flux_linkage_c_list',
        #     shape=(int(s/3),)
        # )
        flux_linkage_a =  self.register_output(
            name='flux_linkage_a',
            var=csdl.sum(delta_A_z_a)
        )
        flux_linkage_b =  self.register_output(
            name='flux_linkage_b',
            var=csdl.sum(delta_A_z_b)
        )
        flux_linkage_c =  self.register_output(
            name='flux_linkage_c',
            var=csdl.sum(delta_A_z_c)
        )

        # flux_linkage_abc =  self.create_output(
        #     name='flux_linkage_abc',
        #     shape=(3,)
        # )
        # flux_linkage_abc[0] = flux_linkage_a
        # flux_linkage_abc[1] = flux_linkage_b[0]
        # flux_linkage_abc[2] = flux_linkage_c[0]
            
        
                
if __name__ == '__main__':
    aaa = FluxLinkageModel()
    sim = Simulator(aaa)

    sim.run()

    print(sim['winding_delta_A_z'])
    print(sim['delta_A_z_a'], sim['delta_A_z_a'].shape)
    print(sim['delta_A_z_b'])
    print(sim['delta_A_z_c'])
    print(sim['flux_linkage_a'].shape)

    
