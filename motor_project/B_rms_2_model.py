# Import `dolfin` first to avoid segmentation fault
from dolfin import *

from csdl import Model, CustomExplicitOperation
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *

class Brms2Model(Model):

    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        A_z = self.declare_variable('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        e = Brms2(fea=self.fea)
        B_rms_2 = csdl.custom(A_z, op=e)
        self.register_output('B_rms_2', B_rms_2)


class Brms2(CustomExplicitOperation):
    """
    input: A_z
    output: B_rms ** 2
    """
    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        self.add_input('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        self.add_output('B_rms_2')
        self.declare_derivatives('*', '*')
        
    def compute(self, inputs, outputs):
        update(self.fea.A_z, inputs['A_z'])
        B_rms_2 = self.fea.extractFluxDensity()
        outputs['B_rms_2'] = B_rms_2

    def compute_derivatives(self, inputs, derivatives):
        update(self.fea.A_z, inputs['A_z'])
        dB_rms_2_dA_z = self.fea.extractFluxDensityDerivatives()
        derivatives['B_rms_2', 'A_z'] = dB_rms_2_dA_z
        
if __name__ == "__main__":
    iq                  = 282.2 / 3
    i_abc               = [
        -iq * np.sin(0.),
        -iq * np.sin(-2*np.pi/3),
        -iq * np.sin(2*np.pi/3),
    ]
    f = open('edge_deformation_data/init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('edge_deformation_data/edge_coord_deltas.txt', 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()
    
    fea = MotorFEA(mesh_file="mesh_files/motor_mesh_1", i_abc=i_abc, 
                            old_edge_coords=old_edge_coords)
    fea.edge_deltas = 0.1*edge_deltas
    sim = Simulator(Brms2Model(fea=fea))
    from matplotlib import pyplot as plt
    print("CSDL: Running the model...")
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim.run()
    print(" B_rms =", np.sqrt(sim['B_rms_2']))
    fea.solveMeshMotion()   
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim.run()
    print(" B_rms =", np.sqrt(sim['B_rms_2']))
    
    print("CSDL: Running check_partials()...")
    sim.check_partials()


        
        
        
        
        
        
        
        
        
        
        
        
        
