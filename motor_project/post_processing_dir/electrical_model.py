import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
# from d_q_transform_model import DQTransformModel

class ElectricalModel(Model): # could also call Battery Model?
    def initialize(self):
        self.parameters.declare('theta')
        self.parameters.declare('flux_linkage_abc')
        self.parameters.declare('flux_linkage_sign')
        self.parameters.declare('instance')

    def define(self):
        theta               = self.parameters['theta']
        flux_linkage_sign   = self.parameters['flux_linkage_sign']  # formatted as [B A C]
        instance            = self.parameters['instance']

        # -------------------- INPUTS TO MODEL --------------------
        p = 12
        # theta           = self.create_input(name='theta', val=0)
        theta           = self.declare_variable(name='theta_{}'.format(instance), val=theta)
        frequency       = self.declare_variable(name='frequency')

        motor_length    = self.declare_variable(name='motor_length', shape=(1,))
        wire_resistance = self.declare_variable('wire_resistance')
        num_windings    = self.declare_variable('num_windings')

        # flux_linkage_abc  = self.create_input(name='flux_linkage_abc', shape=(3,))
        flux_linkage_a_i  = self.declare_variable(name='flux_linkage_a_i_{}'.format(instance), shape=(1,))
        flux_linkage_b_i  = self.declare_variable(name='flux_linkage_b_i_{}'.format(instance), shape=(1,))
        flux_linkage_c_i  = self.declare_variable(name='flux_linkage_c_i_{}'.format(instance), shape=(1,))

        flux_linkage_abc = self.create_output(
            name='flux_linkage_abc_{}'.format(instance),
            shape=(3,)
        )

        flux_linkage_abc[0] = flux_linkage_a_i * motor_length / 2 * num_windings
        flux_linkage_abc[1] = flux_linkage_b_i * motor_length / 2 * num_windings
        flux_linkage_abc[2] = flux_linkage_c_i * motor_length / 2 * num_windings

        # flux_linkage_abc[0] = flux_linkage_a_i * motor_length / 2
        # flux_linkage_abc[1] = flux_linkage_b_i * motor_length / 2
        # flux_linkage_abc[2] = flux_linkage_c_i * motor_length / 2
        
        # -------------------- DQ TRANSFORM MATRICES --------------------
        transform_matrix_forward = self.create_output(
            name='transform_matrix_forward_{}'.format(instance),
            shape=(2,3),
        ) # MOVE THIS VARIABLE OUTSIDE BECAUSE THIS WON'T CHANGE WITH EACH INSTANCE

        transform_matrix_forward[0] = 2/3 * csdl.cos(theta)
        transform_matrix_forward[1] = 2/3 * csdl.cos(theta - 2*np.pi/3)
        transform_matrix_forward[2] = 2/3 * csdl.cos(theta + 2*np.pi/3)
        transform_matrix_forward[3] = -2/3 * csdl.sin(theta)
        transform_matrix_forward[4] = -2/3 * csdl.sin(theta - 2*np.pi/3)
        transform_matrix_forward[5] = -2/3 * csdl.sin(theta + 2*np.pi/3)

        transform_matrix_reverse = self.register_output(
            name='transform_matrix_reverse_{}'.format(instance),
            var= 3/2  *  csdl.transpose(transform_matrix_forward)
        ) # MOVE THIS VARIABLE OUTSIDE BECAUSE THIS WON'T CHANGE WITH EACH INSTANCE

        # -------------------- DQ FLUX LINKAGE --------------------
        flux_linkage_dq = self.register_output(
            name='flux_linkage_dq_{}'.format(instance),
            var=csdl.matvec(transform_matrix_forward, flux_linkage_abc)
        )

        # -------------------- DQ PHASE CURRENT --------------------
        phase_current_dq = self.declare_variable(
            name='phase_current_dq',
            shape=(2,),
        )

        phase_current_abc  = self.register_output(
            name='phase_current_abc_{}'.format(instance), 
            var=csdl.matvec(transform_matrix_reverse, phase_current_dq)
        )

        # -------------------- DQ PHASE VOLTAGE MATRIX --------------------
        phase_voltage_dq = self.create_output(
            name='phase_voltage_dq_{}'.format(instance),
            val=0.,
            shape=(2,)
        )
        phase_voltage_dq[0] = wire_resistance * phase_current_dq[0] - 2 * np.pi * frequency * flux_linkage_dq[1]
        phase_voltage_dq[1] = wire_resistance * phase_current_dq[1] + 2 * np.pi * frequency * flux_linkage_dq[0]

        # -------------------- ABC PHASE VOLTAGE MATRIX --------------------
        phase_voltage_abc = self.create_output(
            name='phase_voltage_abc_{}'.format(instance),
            val=0.,
            shape=(3,)
        )
        phase_voltage_abc[:] = csdl.matvec(transform_matrix_reverse, phase_voltage_dq)

        input_power = self.register_output(
            name='input_power_{}'.format(instance),
            # var=csdl.dot(phase_voltage_abc, phase_current_abc)
            var=3/2 * csdl.dot(phase_voltage_dq, phase_current_dq)
        )

        electromagnetic_torque = self.register_output(
            name='electromagnetic_torque_{}'.format(instance),
            # var=3/2*p/2*(flux_linkage_dq[0]*phase_current_dq[1] - flux_linkage_dq[1]*phase_current_dq[0])
            var=3/2*p/2*(flux_linkage_dq[0]*phase_current_dq[1] - flux_linkage_dq[1]*phase_current_dq[0])
        )



if __name__ == '__main__':
    aaa = ElectricalModel(
        theta = 0,
        flux_linkage_abc = np.array([1., 1., 1.])
    )

    sim = Simulator(aaa)
    sim.run()
    print(sim['transform_matrix_forward'])
    print(sim['transform_matrix_reverse'])
    print(sim['flux_linkage_dq'])
    print(sim['phase_current_dq'])
    print(sim['phase_voltage_dq'])
    print(sim['phase_voltage_abc'])
    print(sim['input_power'])