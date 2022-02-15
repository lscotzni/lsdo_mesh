# Import `dolfin` first to avoid segmentation fault
from dolfin import *

from csdl import Model, CustomExplicitOperation
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *

class AreaModel(Model):

    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size = self.fea.total_dofs_uhat
        uhat = self.declare_variable('uhat',
                        shape=(self.input_size,),
                        val=np.zeros(self.input_size).reshape(self.input_size,))

        e = Area(fea=self.fea)
        winding_area, magnet_area, steel_area = csdl.custom(uhat, op=e)
        self.register_output('winding_area', winding_area)
        self.register_output('magnet_area', magnet_area)
        self.register_output('steel_area', steel_area)


class Area(CustomExplicitOperation):
    """
    input: uhat
    output: winding_area, magnet_area, steel_area
    """
    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size = self.fea.total_dofs_uhat
        self.add_input('uhat',
                        shape=(self.input_size,),
                        val=0.0)
        self.add_output('winding_area')
        self.add_output('magnet_area')
        self.add_output('steel_area')
        self.declare_derivatives('*', '*')
        
    def compute(self, inputs, outputs):
        update(self.fea.uhat, inputs['uhat'])
        winding_area, magnet_area, steel_area = self.fea.calcAreas()
        outputs['winding_area'] = winding_area
        outputs['magnet_area'] = magnet_area
        outputs['steel_area'] = steel_area

    def compute_derivatives(self, inputs, derivatives):
        update(self.fea.uhat, inputs['uhat'])
        dWA, dMA, dSA = self.fea.calcAreaDerivatives()
        derivatives['winding_area', 'uhat'] = dWA
        derivatives['magnet_area', 'uhat'] = dMA
        derivatives['steel_area', 'uhat'] = dSA
        
if __name__ == "__main__":
    iq                  = 282.2 / 3
    f = open('coarse_mesh_Ru/init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('coarse_mesh_Ru/edge_coord_deltas.txt', 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()
    
#    fea = MotorFEA(mesh_file="mesh_files/motor_mesh_1", 
    fea = MotorFEA(mesh_file="coarse_mesh_Ru/motor_mesh_coarse_1",
                            old_edge_coords=old_edge_coords)
    fea.edge_deltas = 0.1*edge_deltas
    updateR(fea.iq, iq)
    sim = Simulator(AreaModel(fea=fea))
    
    from matplotlib import pyplot as plt
    print("CSDL: Running the model...")
    
    fea.solveMeshMotion()   
    sim['uhat'] = fea.uhat.vector().get_local()
    sim.run()
    print(sim['winding_area'])
    print(sim['magnet_area'])
    print(sim['steel_area'])
    
    print("CSDL: Running check_partials()...")
    sim.check_partials()


        
        
        
        
        
        
        
        
        
        
        
        
        
