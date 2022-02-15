# Import `dolfin` first to avoid segmentation fault
from dolfin import *

from csdl import Model, CustomExplicitOperation
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *

class FluxInfluencePMModel(Model):

    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size_A_z = self.fea.total_dofs_A_z
        self.input_size_uhat = self.fea.total_dofs_uhat
        A_z = self.declare_variable('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        uhat = self.declare_variable('uhat',
                        shape=(self.input_size_uhat,),
                        val=0.0)
        e = FluxInfluencePM(fea=self.fea)
        B_influence_pm = csdl.custom(A_z, uhat, op=e)
        self.register_output('B_influence_pm', B_influence_pm)


class FluxInfluencePM(CustomExplicitOperation):
    """
    input: A_z, uhat
    output: B_influence_pm = B**2*dx(subdomains)
    """
    def initialize(self):
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        # parameters for Eddy current loss
        self.subdomains = self.fea.pm_loss_subdomain
        self.exponent = 2
        
        self.input_size_A_z = self.fea.total_dofs_A_z
        self.input_size_uhat = self.fea.total_dofs_uhat
        A_z = self.add_input('A_z',
                        shape=(self.input_size_A_z,),
                        val=0.0)
        uhat = self.add_input('uhat',
                        shape=(self.input_size_uhat,),
                        val=0.0)
        self.add_output('B_influence_pm')
        self.declare_derivatives('*', '*')
        
    def compute(self, inputs, outputs):
        update(self.fea.A_z, inputs['A_z'])
        update(self.fea.uhat, inputs['uhat'])
        B_influence_pm = self.fea.calcFluxInfluence(
                        n=self.exponent, subdomains=self.subdomains)
        outputs['B_influence_pm'] = B_influence_pm

    def compute_derivatives(self, inputs, derivatives):
        update(self.fea.A_z, inputs['A_z'])
        update(self.fea.uhat, inputs['uhat'])
        dFdAz, dFduhat = self.fea.calcFluxInfluenceDerivatives(
                        n=self.exponent, subdomains=self.subdomains)
        derivatives['B_influence_pm', 'A_z'] = dFdAz
        derivatives['B_influence_pm', 'uhat'] = dFduhat
        
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
    
    fea = MotorFEA(mesh_file="mesh_files/motor_mesh_1",
                            old_edge_coords=old_edge_coords)
    fea.edge_deltas = 0.1*edge_deltas
    fea.iq.assign(Constant(float(iq)))
    sim = Simulator(FluxInfluencePMModel(fea=fea))
    from matplotlib import pyplot as plt
    print("CSDL: Running the model...")
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim['uhat'] = fea.uhat.vector().get_local()
    sim.run()
    print(" B_influence_pm =", sim['B_influence_pm'])
    fea.solveMeshMotion()   
    fea.solveMagnetostatic()
    sim['A_z'] = fea.A_z.vector().get_local()
    sim['uhat'] = fea.uhat.vector().get_local()
    sim.run()
    print(" B_influence_pm =", sim['B_influence_pm'])
    
    print("CSDL: Running check_partials()...")
    sim.check_partials()


        
        
        
        
        
        
        
        
        
        
        
        
        
