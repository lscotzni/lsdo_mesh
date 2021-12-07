import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator
from d_q_transform_model import DQTransformModel

class FluxLinkageModel(Model): # could also call Battery Model?
    def initialize(self):
        self.parameters.declare('theta')
        self.parameters.declare('flux_linkage_abc')

    def define(self):
        # theta               = self.parameters['theta']
        # flux_linkage_abc    = self.parameters['flux_linkage_abc']

        theta  = self.create_input(
            name='theta',
            val=0.,

        )

        flux_linkage_abc  = self.create_input(
            name='flux_linkage_abc',
            shape=(3,)
        )

        phase_current_abc  = self.create_input(
            name='phase_current_abc',
            shape=(3,)
        )

        wire_resistance = self.create_input(
            name='wire_resistance',
            shape=(1,1)
        )

        # transform_matrix_forward = 2/3 * np.array([
        #     np.cos(theta), 
        #     np.cos(theta - 2*np.pi/3), 
        #     np.cos(theta + 2*np.pi/3),
        #     -np.sin(theta),
        #     -np.sin(theta - 2*np.pi/3),
        #     -np.sin(theta + 2*np.pi/3),
        # ]).reshape((2,3))

        # transform_matrix_forward = self.declare_variable(
        #     name='transform_matrix_forward',
        #     val=transform_matrix_forward,
        #     shape=transform_matrix_forward.shape
        # )

        transform_matrix_forward = self.create_output(
            name='transform_matrix_forward',
            shape=(2,3)
        )

        transform_matrix_forward[0] = csdl.cos(theta)
        transform_matrix_forward[1] = csdl.cos(theta - 2*np.pi/3)
        transform_matrix_forward[2] = csdl.cos(theta + 2*np.pi/3)
        transform_matrix_forward[3] = -csdl.sin(theta)
        transform_matrix_forward[4] = -csdl.sin(theta - 2*np.pi/3)
        transform_matrix_forward[5] = -csdl.sin(theta + 2*np.pi/3)
        
        
        # transform_matrix_reverse = np.array([
        #     np.cos(theta), 
        #     -np.sin(theta),
        #     np.cos(theta - 2*np.pi/3),
        #     -np.sin(theta - 2*np.pi/3), 
        #     np.cos(theta + 2*np.pi/3),
        #     -np.sin(theta + 2*np.pi/3),
        # ]).reshape((3,2))
        # transform_matrix_reverse = self.declare_variable(
        #     name='transform_matrix_reverse',
        #     val=transform_matrix_reverse,
        #     shape=transform_matrix_reverse.shape
        # )

        # flux_linkage_abc = self.declare_variable(
        #     name='flux_linkage_abc',
        #     val=flux_linkage_abc,
        #     shape=flux_linkage_abc.shape
        # )

        flux_linkage_dq = self.register_output(
            name='flux_linkage_dq',
            var=csdl.matvec(transform_matrix_forward, flux_linkage_abc)
        )

        phase_current_dq = self.register_output(
            name='phase_current_dq',
            var=csdl.matvec(transform_matrix_forward, phase_current_abc)
        )

        # resistance_matrix_dq = self.create_output(
        #     name='resistance_matrix_dq',
        #     shape=(2,2)
        # )

        # resistance_matrix_dq[0,1] = wire_resistance
        # resistance_matrix_dq[1,0] = wire_resistance


if __name__ == '__main__':
    aaa = FluxLinkageModel(
        theta = 0,
        flux_linkage_abc = np.array([1., 1., 1.])
    )

    sim = Simulator(aaa)
    print(sim['transform_matrix_forward'])
    print(sim['flux_linkage_dq'])
    print(sim['phase_current_dq'])
    # print(sim['resistance_matrix_dq'])