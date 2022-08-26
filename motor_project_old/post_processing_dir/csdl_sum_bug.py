from csdl_om import Simulator
from csdl import Model
import csdl
import numpy as np


class ExampleSingleVector(Model):
    def define(self):
        n = 5

        # Declare a vector of length 3 as input
        v1 = self.declare_variable('v1', val=np.arange(n))
        v2 = self.declare_variable('v2', val=np.arange(n-1))

        # Output the sum of all the elements of the vector v1
        single_vector_sum_1 =  self.register_output('single_vector_sum_1', csdl.sum(v1, axes=(0,)))
        single_vector_sum_2 =  self.register_output('single_vector_sum_2', csdl.sum(v2, axes=(0,)))

        sum_vector = self.create_output(
            name='sum_vector',
            shape=(2,)
        )

        print(sum_vector.shape)
        print(single_vector_sum_1.shape)
        print(single_vector_sum_2.shape)

        sum_vector[0] = single_vector_sum_1
        sum_vector[1] = single_vector_sum_2


sim = Simulator(ExampleSingleVector())
sim.run()

print('v1', sim['v1'].shape)
print(sim['v1'])
print('single_vector_sum_1', sim['single_vector_sum_1'].shape)
print(sim['single_vector_sum_2'])
print('single_vector_sum_2', sim['single_vector_sum_2'].shape)
print(sim['single_vector_sum_2'])
print('sum_vector', sim['sum_vector'].shape)
print(sim['sum_vector'])