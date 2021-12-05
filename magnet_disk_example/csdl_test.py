from csdl import Model, CustomImplicitOperation, NewtonSolver, ScipyKrylov
import csdl
import numpy as np

class M(Model):
    def define(self):
        a = self.declare_variable('a', val=1.)
        b = self.declare_variable('b', val=-4.)
        c = self.declare_variable('c', val=3.)

#        x = csdl.custom(a, b, c, op=ExampleImplicitSimple())

        e = ExampleImplicitSimple()
        e.linear_solver = ScipyKrylov()
        e.nonlinear_solver = NewtonSolver(solve_subsystems=False)
        x = csdl.custom(a, b, c, op=e)
        self.register_output('x', x)

class ExampleImplicitSimple(CustomImplicitOperation):
    """
    :param var: x
    """
    def define(self):
        print("="*40)
        print(" Running define()...")
        print("="*40)
        self.add_input('a', val=1.)
        self.add_input('b', val=-4.)
        self.add_input('c', val=3.)
        self.add_output('x', val=0.)
        self.declare_derivatives('x', 'x')
        self.declare_derivatives('x', ['a', 'b', 'c'])

#        self.linear_solver = ScipyKrylov()
#        self.nonlinear_solver = NewtonSolver(solve_subsystems=False)

    def evaluate_residuals(self, inputs, outputs, residuals):
        print("="*40)
        print(" Running evaluate_residual()...")
        print("="*40)
        x = outputs['x']
        a = inputs['a']
        b = inputs['b']
        c = inputs['c']
        residuals['x'] = a * x**2 + b * x + c

    def compute_derivatives(self, inputs, outputs, derivatives):
        print("="*40)
        print(" Running compute_derivatives()...")
        print("="*40)
        a = inputs['a']
        b = inputs['b']
        x = outputs['x']

        derivatives['x', 'a'] = x**2
        derivatives['x', 'b'] = x
        derivatives['x', 'c'] = 1.0
        derivatives['x', 'x'] = 2 * a * x + b
        
from csdl_om import Simulator

# Generate an implementation.
sim = Simulator(M())
# Run simulation.
#sim.run()
sim.check_partials()
# Access values
print(sim['x'])
