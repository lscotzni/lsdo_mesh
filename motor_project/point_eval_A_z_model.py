# Import `dolfin` first to avoid segmentation fault
from dolfin import *

from csdl import Model, CustomExplicitOperation
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *

class PointEvalAzModel(Model):

    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        A_z = self.declare_variable('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        e = PointEvalAz(fea=self.fea)
        A_p = csdl.custom(A_z, op=e)
        self.register_output('A_p', A_p)


class PointEvalAz(CustomExplicitOperation):
    """
    input: A_z
    output: point evaluations of A_z
    """
    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.points = self.fea.getPointCoords()
        self.input_size_A_z = self.fea.total_dofs_A_z
        self.add_input('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        self.add_output('A_p',
                        shape=(36,),
                        val=0.0)
        self.declare_derivatives('A_p', 'A_z')
        
    def compute(self, inputs, outputs):
        update(self.fea.A_z, inputs['A_z'])
        outputs['A_p'] = self.fea.evalPointAz(self.points)

    def compute_derivatives(self, inputs, derivatives):
        update(self.fea.A_z, inputs['A_z'])
        dApdAz = self.fea.getPointEvalDerivatives(self.points)
        derivatives['A_p', 'A_z'] = dApdAz.todense()

                    
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
    fea.edge_deltas = edge_deltas
    sim = Simulator(PointEvalAzModel(fea=fea))
    from matplotlib import pyplot as plt
    print("CSDL: Running the model...")
    
    
#    fea.solveMeshMotion()   
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim.run()
    print("- "*30)
    print("Point evaluations of A_z along a radius of 0.05:")
    print("- "*30)
    print(sim['A_p'])
    
#    print("CSDL: Running check_partials()...")
#    sim.check_partials()


        
        
        
        
        
        
        
        
        
        
        
        
        
