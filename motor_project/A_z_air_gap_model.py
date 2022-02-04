# Import `dolfin` first to avoid segmentation fault
from dolfin import *

from csdl import Model, CustomExplicitOperation
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *

class AzAirGapModel(Model):

    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        A_z = self.declare_variable('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        e = AzAirGap(fea=self.fea)
        A_z_air_gap = csdl.custom(A_z, op=e)
        self.register_output('A_z_air_gap', A_z_air_gap)


class AzAirGap(CustomExplicitOperation):
    """
    input: A_z
    output: nodal evaluations of A_z at the air gap
    """
    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        self.output_size = len(self.fea.A_z_air_gap_indices)
        self.add_input('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        self.add_output('A_z_air_gap',
                        shape=(self.output_size,),
                        val=0.0)
        self.declare_derivatives('A_z_air_gap', 'A_z')
        
    def compute(self, inputs, outputs):
        update(self.fea.A_z, inputs['A_z'])
        outputs['A_z_air_gap'] = self.fea.extractAzAirGap()

    def compute_derivatives(self, inputs, derivatives):
        update(self.fea.A_z, inputs['A_z'])
        dA_ag_dA_z = self.fea.extractAzAirGapDerivatives()
        derivatives['A_z_air_gap', 'A_z'] = dA_ag_dA_z.todense()

                    
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
    fea.A_z_air_gap_indices = np.array([1,3,5,7,9])
    
    sim = Simulator(AzAirGapModel(fea=fea))
    from matplotlib import pyplot as plt
    print("CSDL: Running the model...")
#    fea.solveMeshMotion()   
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim.run()
    print("- "*30)
    print("Nodal evaluations of A_z in the air gap:")
    print("- "*30)
    print(sim['A_z_air_gap'])
    
    print("CSDL: Running check_partials()...")
    sim.check_partials()


        
        
        
        
        
        
        
        
        
        
        
        
        
