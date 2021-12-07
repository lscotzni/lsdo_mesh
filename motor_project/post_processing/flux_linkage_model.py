import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class FluxLinkageModel(Model):
    def define(self):
        s = 36
        motor_length = self.create_input(name='motor_length')
        num_windings = self.create_input(name='num_windings')

        mag_vec_pot_stator_teeth = self.create_input(
            name='mag_vec_pot_stator_teeth',
            shape=(2*s,)
        )
        # THIS HOLDS THE AVERAGE MAGNETIC VECTOR POTENTIAL VALUE IN A WINDING
        # WE NEED TO TAKE THE DELTA

        mag_vec_pot_a_list =  self.create_input(
            name='mag_vec_pot_a_list',
            shape=(int(2*s/3),)
        )

        mag_vec_pot_b_list =  self.create_input(
            name='mag_vec_pot_b_list',
            shape=(int(2*s/3),)
        )

        mag_vec_pot_c_list =  self.create_input(
            name='mag_vec_pot_c_list',
            shape=(int(2*s/3),)
        )


        flux_linkage_a_list =  self.create_output(
            name='flux_linkage_a_list',
            shape=(int(s/3),)
        )

        flux_linkage_b_list =  self.create_output(
            name='flux_linkage_b_list',
            shape=(int(s/3),)
        )

        flux_linkage_c_list =  self.create_output(
            name='flux_linkage_c_list',
            shape=(int(s/3),)
        )

        for i in range(int(s/3)): # need absolute values around these; might be better to find absolute deltas outside
            flux_linkage_a_list[i] = mag_vec_pot_a_list[2*i + 1] - mag_vec_pot_a_list[2*i]
            flux_linkage_b_list[i] = mag_vec_pot_b_list[2*i + 1] - mag_vec_pot_b_list[2*i]
            flux_linkage_c_list[i] = mag_vec_pot_c_list[2*i + 1] - mag_vec_pot_c_list[2*i]

        if False:
            for i in range(int(s/3)):
                if (i+1)%3 == 1: # PHASE A
                    for j in range(3):
                        flux_linkage_a_list[i] = 1
                    # 6*i to 6*i + 6
                    
                elif (i+1)%3 == 2: # PHASE B
                    for j in range(3):
                        flux_linkage_a_list[i] = 1
                    # 6*i to 6*i + 6
                    
                elif (i+1)%3 == 0: # PHASE C
                    for j in range(3):
                        flux_linkage_a_list[i] = 1
                    # 6*i to 6*i + 6
                


