import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

class DQTransformModel(Model):
    def initialize(self):
        self.parameters.declare('theta')
        self.parameters.declare('direction')
        self.parameters.declare('input_vector')
    def define(self):
        theta = self.parameters['theta']
        direction = self.parameters['direction']
        input_vector = self.parameters['input_vector']

        if direction == 'forward': # abc -> dq
            transform_matrix = 2/3 * np.array([
                np.cos(theta), 
                np.cos(theta - 2*np.pi/3), 
                np.cos(theta + 2*np.pi/3),
                -np.sin(theta),
                -np.sin(theta - 2*np.pi/3),
                -np.sin(theta + 2*np.pi/3),
            ]).reshape((2,3))

        elif direction == 'reverse': # dq -> abc
            transform_matrix = np.array([
                np.cos(theta), 
                -np.sin(theta),
                np.cos(theta - 2*np.pi/3),
                -np.sin(theta - 2*np.pi/3), 
                np.cos(theta + 2*np.pi/3),
                -np.sin(theta + 2*np.pi/3),
            ]).reshape((3,2))
            
        print(transform_matrix) 

        in_vector = self.declare_variable(
            name='input_vector',
            val=input_vector,
            shape=input_vector.shape
        )

        transform_matrix = self.declare_variable(
            name='transform_matrix',
            val=transform_matrix,
            shape=transform_matrix.shape
        )

        out_vector = self.register_output(
            name='output_vector',
            var=csdl.matvec(transform_matrix, in_vector),
            # shape=(transform_matrix.shape[0],)
        )
        
        # self.register_output(
        #     'output_vector',
        #     var=out_vector
        # )
        
if __name__ == '__main__':
    aaa = DQTransformModel(
        theta=0,
        direction='forward',
        input_vector=np.array([1, 1, 1])
    )
    sim = Simulator(aaa)
    print(sim['output_vector'])

    bbb = DQTransformModel(
        theta=0,
        direction='reverse',
        input_vector=np.array([1, 1])
    )
    sim = Simulator(bbb)
    print(sim['output_vector'])