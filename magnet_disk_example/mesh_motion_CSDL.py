from csdl import CustomImplicitOperation, Model
import csdl
import numpy as np
from csdl_om import Simulator
from magnetostatic_fea import *


class M(Model):

    def initialize(self):
        print("="*40)
        print("CSDL: Running initialize()...")
        print("="*40)
        self.parameters.declare('fea')
        
    def define(self):
        self.fea = self.parameters['fea']
        self.input_size = len(self.fea.edge_indices)
        self.output_size = self.fea.total_dofs_uhat
        edge_deltas = self.declare_variable('edge_deltas', 
                        shape=(self.input_size,), 
                        val=0.1*np.ones(self.input_size).reshape(self.input_size,))

        e = MeshMotion(fea=self.fea)
        uhat = csdl.custom(edge_deltas, op=e)
        self.register_output('uhat', uhat)
        
        
class MeshMotion(CustomImplicitOperation):
    """
    input: edge_deltas
    output: uhat
    """
    def initialize(self):
        print("="*40)
        print("CSDL: Running initialize()...")
        print("="*40)
        self.parameters.declare('fea')
        
    def define(self):
        print("="*40)
        print("CSDL: Running define()...")
        print("="*40)
        self.fea = self.parameters['fea']
        self.input_size = len(self.fea.edge_indices)
        self.output_size = self.fea.total_dofs_uhat
        self.add_input('edge_deltas', 
                        shape=(self.input_size,), 
                        val=0.1*np.ones(self.input_size).reshape(self.input_size,))
        self.add_output('uhat',
                        shape=(self.output_size,),
                        val=0.1*np.ones(self.output_size).reshape(self.output_size,))
        self.declare_derivatives('uhat', 'uhat')       
        self.declare_derivatives('uhat', 'edge_deltas')
     
    def updateMeshMotionBC(self, inputs):
        self.fea.edge_deltas = inputs['edge_deltas']
        self.fea.uhat_old.vector()[self.fea.edge_indices] = inputs['edge_deltas']
        self.bcs = self.fea.setBCMeshMotion()
        
    def evaluate_residuals(self, inputs, outputs, residuals):
        print("="*40)
        print("CSDL: Running evaluate_residuals()...")
        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])
        resM = assemble(self.fea.resM())
        uhat_0 = Function(self.fea.VHAT)
        uhat_1 = Function(self.fea.VHAT)
        uhat_1.assign(self.fea.uhat)
        uhat_1.vector()[self.fea.edge_indices] = inputs['edge_deltas']
        uhat_0.assign(self.fea.uhat-uhat_1)
        
        self.bcs.apply(resM, uhat_0.vector())
        
        residuals['uhat'] = resM.get_local()
    
    def solve_residual_equations(self, inputs, outputs):
        print("="*40)
        print("CSDL: Running solve_residual_equations()...")
        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])
        
        self.fea.solveMeshMotion()
               
        outputs['uhat'] = self.fea.uhat.vector().get_local()
        update(self.fea.uhat, outputs['uhat'])
        
        
    def compute_derivatives(self, inputs, outputs, derivatives):
        print("="*40)
        print("CSDL: Running compute_derivatives()...")
        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])
        
        self.dRdu = assemble(self.fea.dRm_duhat)
        self.bcs.apply(self.dRdu)
        self.dRdf = convertToDense(self.fea.getBCDerivatives())
        self.A = self.dRdu
        for i in range(len(self.fea.edge_indices)):
            row = self.fea.edge_indices[i].astype('int')
            col = i
            print(row, col)
            print(self.dRdf[row][col])
        self.niter = 0
        
    def compute_jacvec_product(self, inputs, outputs, 
                                d_inputs, d_outputs, d_residuals, mode):
        print("="*40)
        print("CSDL: Running compute_jacvec_product()...")
        print("="*40)
#        self.niter += 1
#        if self.niter >= 4:
#            exit()
        if mode == 'fwd':
            if 'uhat' in d_residuals:
                if 'uhat' in d_outputs:
                    update(self.fea.duhat, d_outputs['uhat'])
                    d_residuals['uhat'] += computeMatVecProductFwd(
                                            self.dRdu, self.fea.duhat)
                if 'edge_deltas' in d_inputs:
                    d_residuals['uhat'] += np.dot(
                                            self.dRdf, 
                                            d_inputs['edge_deltas'])
                    
        if mode == 'rev':
            if 'uhat' in d_residuals:
                update(self.fea.dRhat, d_residuals['uhat'])
                if 'uhat' in d_outputs:
                    d_outputs['uhat'] += computeMatVecProductBwd(
                                            self.dRdu, self.fea.dRhat)
                if 'edge_deltas' in d_inputs:
                    d_inputs['edge_deltas'] += np.dot(
                                            np.transpose(self.dRdf), 
                                            d_residuals['uhat'])
                    
    def apply_inverse_jacobian(self, d_outputs, d_residuals, mode):
        print("="*40)
        print("CSDL: Running apply_inverse_jacobian()...")
        print("="*40)
        
        if mode == 'fwd':
            d_outputs['uhat'] = self.fea.solveLinearFwd(
                                                self.A, 
                                                d_residuals['uhat'])
        else:
            d_residuals['uhat'] = self.fea.solveLinearBwd(
                                                self.A, 
                                                d_outputs['uhat'])
            
            
if __name__ == "__main__":
    fea = MagnetostaticProblem()
    sim = Simulator(M(fea=fea))
    
#    sim['edge_deltas'] = generateMeshMovement(-pi/36)
#    sim.run()
    
#    fea.uhat.vector().set_local(sim['uhat'])
#    plt.figure(1)
#    fea.moveMesh()
#    plot(fea.mesh)
#    plt.show()
    print(sim['uhat'].shape, sim['edge_deltas'].shape)
    
    # TODO: the partial derivative of residual wrt edge_deltas is not correct
    sim.check_partials()
