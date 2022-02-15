import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class TimeAverageModel(Model):
    def initialize(self):
        self.parameters.declare('instances')
        
    def define(self):
        instances   = self.parameters['instances']
        flux_influence_ec_list  = self.create_output(
            name='flux_influence_ec_list',
            shape=(instances,)
        )

        flux_influence_h_list  = self.create_output(
            name='flux_influence_h_list',
            shape=(instances,)
        )

        flux_influence_pm_list  = self.create_output(
            name='flux_influence_pm_list',
            shape=(instances,)
        )

        em_torque_list  = self.create_output(
            name='em_torque_list',
            shape=(instances,)
        )
        
        for i in range(instances):
            temp_flux_influence_ec  = self.declare_variable('B_influence_ec_{}'.format(i+1))
            flux_influence_ec_list[i] = temp_flux_influence_ec

            temp_flux_influence_h   = self.declare_variable('B_influence_hysteresis_{}'.format(i+1))
            flux_influence_h_list[i] = temp_flux_influence_h

            temp_flux_influence_pm  = self.declare_variable('B_influence_pm_{}'.format(i+1))
            flux_influence_pm_list[i] = temp_flux_influence_pm

            temp_em_torque = self.declare_variable('electromagnetic_torque_{}'.format(i+1))
            em_torque_list[i] = temp_em_torque

        avg_flux_influence_ec   = csdl.average(flux_influence_ec_list)
        avg_flux_influence_ec   = self.register_output(
            name='avg_flux_influence_ec',
            var=avg_flux_influence_ec
        )

        avg_flux_influence_h    = csdl.average(flux_influence_h_list)
        avg_flux_influence_h   = self.register_output(
            name='avg_flux_influence_h',
            var=avg_flux_influence_h
        )

        avg_flux_influence_pm   = csdl.average(flux_influence_pm_list)
        avg_flux_influence_pm   = self.register_output(
            name='avg_flux_influence_pm',
            var=avg_flux_influence_pm
        )

        avg_electromagnetic_torque  = csdl.average(em_torque_list)
        avg_electromagnetic_torque  = self.register_output(
            name='avg_electromagnetic_torque',
            var=avg_electromagnetic_torque
        )

            
if __name__  == '__main__':
    aaa = TimeAverageModel()
    sim = Simulator(aaa)
    sim.run()
    # sim.visualize_model()