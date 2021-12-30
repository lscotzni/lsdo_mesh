import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
# from d_q_transform_model import DQTransformModel

class ElectricalModel(Model): # could also call Battery Model?
    def initialize(self):
        self.parameters.declare('theta')
        self.parameters.declare('flux_linkage_abc')

    def define(self):
        # theta               = self.parameters['theta']
        # flux_linkage_abc    = self.parameters['flux_linkage_abc']

        # -------------------- INPUTS TO MODEL --------------------
        p = 12
        torque_scale = self.create_input(name='torque_scale', val=3/2*p/2)
        # theta  = self.create_input(name='theta', val=0.)
        theta  = self.create_input(name='theta', val=np.pi/2)
        omega   = self.declare_variable(name='omega')
        # phase_current_abc  = self.create_input(
        #     name='phase_current_abc', 
        #     val=np.array([94.26 * 39, 47.13 * 39, -47.13 * 39]),
        #     # val=np.array([66.65 * 39, -91.04 * 39, 24.4 * 39]),
        #     shape=(3,)
        # )

        motor_length = self.declare_variable(name='motor_length', shape=(1,))
        num_windings = self.declare_variable(name='num_windings', shape=(1,), val=39)
        winding_area = self.declare_variable(name='winding_area')

        # flux_linkage_abc  = self.create_input(name='flux_linkage_abc', shape=(3,))
        flux_linkage_a_i  = self.declare_variable(name='flux_linkage_a_i', shape=(1,))
        flux_linkage_b_i  = self.declare_variable(name='flux_linkage_b_i', shape=(1,))
        flux_linkage_c_i  = self.declare_variable(name='flux_linkage_c_i', shape=(1,))

        flux_linkage_abc = self.create_output(
            name='flux_linkage_abc',
            shape=(3,)
        )
        flux_linkage_abc[0] = flux_linkage_a_i * motor_length * num_windings
        flux_linkage_abc[1] = flux_linkage_b_i * motor_length * num_windings
        flux_linkage_abc[2] = flux_linkage_c_i * motor_length * num_windings
        
        # wire_resistance = self.create_input(name='wire_resistance')
        # wire_resistance = 1.68e-8 * motor_length / (winding_area/39) * 72
        # wire_resistance = 1.68e-8 * motor_length / (winding_area/39) * 18
        wire_resistance = 1.68e-8 * motor_length / (winding_area/num_windings) * 18 * num_windings
        wire_resistance = self.register_output(
            name='wire_resistance',
            var = wire_resistance
        )

        # -------------------- DQ TRANSFORM MATRICES --------------------
        transform_matrix_forward = self.create_output(
            name='transform_matrix_forward',
            shape=(2,3),
        )

        transform_matrix_forward[0] = 2/3 * csdl.cos(theta)
        transform_matrix_forward[1] = 2/3 * csdl.cos(theta - 2*np.pi/3)
        transform_matrix_forward[2] = 2/3 * csdl.cos(theta + 2*np.pi/3)
        transform_matrix_forward[3] = -2/3 * csdl.sin(theta)
        transform_matrix_forward[4] = -2/3 * csdl.sin(theta - 2*np.pi/3)
        transform_matrix_forward[5] = -2/3 * csdl.sin(theta + 2*np.pi/3)

        transform_matrix_reverse = self.register_output(
            name='transform_matrix_reverse',
            var= 3/2  *  csdl.transpose(transform_matrix_forward)
        )

        # -------------------- DQ FLUX LINKAGE --------------------
        flux_linkage_dq = self.register_output(
            name='flux_linkage_dq',
            var=csdl.matvec(transform_matrix_forward, flux_linkage_abc)
        )

        # -------------------- DQ PHASE CURRENT --------------------
        # phase_current_dq = self.register_output(
        #     name='phase_current_dq',
        #     var=csdl.matvec(transform_matrix_forward, phase_current_abc)
        # )
        phase_current_dq = self.declare_variable(
            name='phase_current_dq',
            shape=(2,),
        )

        phase_current_abc  = self.register_output(
            name='phase_current_abc', 
            var=csdl.matvec(transform_matrix_reverse, phase_current_dq)
        )

        # -------------------- DQ PHASE VOLTAGE MATRIX --------------------
        phase_voltage_dq = self.create_output(
            name='phase_voltage_dq',
            val=0.,
            shape=(2,)
        )
        phase_voltage_dq[0] = wire_resistance * phase_current_dq[0] - omega * flux_linkage_dq[1]
        phase_voltage_dq[1] = wire_resistance * phase_current_dq[1] + omega * flux_linkage_dq[0]

        # -------------------- ABC PHASE VOLTAGE MATRIX --------------------
        phase_voltage_abc = self.create_output(
            name='phase_voltage_abc',
            val=0.,
            shape=(3,)
        )
        phase_voltage_abc[:] = csdl.matvec(transform_matrix_reverse, phase_voltage_dq)

        input_power = self.register_output(
            name='input_power',
            # var=csdl.dot(phase_voltage_abc, phase_current_abc)
            var=3/2 * csdl.dot(phase_voltage_dq, phase_current_dq)
        )

        output_torque = self.register_output(
            name='output_torque',
            var=torque_scale*(flux_linkage_dq[0]*phase_current_dq[1] - flux_linkage_dq[1]*phase_current_dq[0])
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