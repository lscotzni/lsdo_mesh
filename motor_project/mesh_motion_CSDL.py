from dolfin import *
from csdl import CustomImplicitOperation, Model
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *


class M(Model):

    def initialize(self):
#        print("="*40)
#        print("CSDL: Running initialize()...")
#        print("="*40)
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size = len(self.fea.edge_indices)
        self.output_size = self.fea.total_dofs_uhat
        edge_deltas = self.declare_variable('edge_deltas',
                        shape=(self.input_size,),
                        val=0.1*np.ones(self.input_size).reshape(self.input_size,)
                        )

        e = MeshMotion(fea=self.fea)
        uhat = csdl.custom(edge_deltas, op=e)
        self.register_output('uhat', uhat)


class MeshMotion(CustomImplicitOperation):
    """
    input: edge_deltas
    output: uhat
    """
    def initialize(self):
#        print("="*40)
#        print("CSDL: Running initialize()...")
#        print("="*40)
        self.parameters.declare('fea')

    def define(self):
#        print("="*40)
#        print("CSDL: Running define()...")
#        print("="*40)
        self.fea = self.parameters['fea']
        self.input_size = len(self.fea.edge_indices)
        self.output_size = self.fea.total_dofs_uhat
        self.add_input('edge_deltas',
                        shape=(self.input_size,),
                        val=0.1*np.ones(self.input_size).reshape(self.input_size,)
                        )
        self.add_output('uhat',
                        shape=(self.output_size,),
                        val=0.1*np.ones(self.output_size).reshape(self.output_size,)
                        )
        self.declare_derivatives('uhat', 'uhat')
        self.declare_derivatives('uhat', 'edge_deltas')

    def updateMeshMotionBC(self, inputs):
        self.fea.edge_deltas = inputs['edge_deltas']
        for i in range(self.input_size):
            self.fea.uhat_old.vector()[self.fea.edge_indices[i]] = inputs['edge_deltas'][i]
        self.bcs = self.fea.setBCMeshMotion()

    def evaluate_residuals(self, inputs, outputs, residuals):
#        print("="*40)
#        print("CSDL: Running evaluate_residuals()...")
#        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])
        resM = assemble(self.fea.resM())
        uhat_0 = Function(self.fea.VHAT)
        uhat_0.assign(self.fea.uhat)

        self.bcs.apply(resM, uhat_0.vector())
        residuals['uhat'] = resM.get_local()

    def solve_residual_equations(self, inputs, outputs):
#        print("="*40)
#        print("CSDL: Running solve_residual_equations()...")
#        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])

        self.fea.solveMeshMotion()

        outputs['uhat'] = self.fea.uhat.vector().get_local()
        update(self.fea.uhat, outputs['uhat'])


    def compute_derivatives(self, inputs, outputs, derivatives):
#        print("="*40)
#        print("CSDL: Running compute_derivatives()...")
#        print("="*40)
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])

        self.dRdu = assemble(self.fea.dRm_duhat)
        self.bcs.apply(self.dRdu)
        self.dRdf = self.fea.getBCDerivatives()
        self.A = self.dRdu


    def compute_jacvec_product(self, inputs, outputs,
                                d_inputs, d_outputs, d_residuals, mode):
#        print("="*40)
#        print("CSDL: Running compute_jacvec_product()...")
#        print("="*40)
        if mode == 'fwd':
            if 'uhat' in d_residuals:
                if 'uhat' in d_outputs:
                    update(self.fea.duhat, d_outputs['uhat'])
                    d_residuals['uhat'] += computeMatVecProductFwd(
                                            self.dRdu, self.fea.duhat)
                if 'edge_deltas' in d_inputs:
                    d_residuals['uhat'] += self.dRdf.dot(d_inputs['edge_deltas'])

        if mode == 'rev':
            if 'uhat' in d_residuals:
                update(self.fea.dRhat, d_residuals['uhat'])
                if 'uhat' in d_outputs:
                    d_outputs['uhat'] += computeMatVecProductBwd(
                                            self.dRdu, self.fea.dRhat)
                if 'edge_deltas' in d_inputs:
                    d_inputs['edge_deltas'] += self.dRdf.transpose().dot(
                                            d_residuals['uhat'])

    def apply_inverse_jacobian(self, d_outputs, d_residuals, mode):
#        print("="*40)
#        print("CSDL: Running apply_inverse_jacobian()...")
#        print("="*40)

        if mode == 'fwd':
            d_outputs['uhat'] = self.fea.solveLinearFwd(
                                                self.A,
                                                d_residuals['uhat'])
        else:
            d_residuals['uhat'] = self.fea.solveLinearBwd(
                                                self.A,
                                                d_outputs['uhat'])


if __name__ == "__main__":
    iq                  = 282.2 / 3
    i_abc               = [
        -iq * np.sin(0.),
        -iq * np.sin(-2*np.pi/3),
        -iq * np.sin(2*np.pi/3),
    ]
    f = open('init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('edge_coord_deltas.txt', 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()
    # One-time computation for the initial edge coordinates from
    # the code that creates the mesh file
    # old_edge_coords = getInitialEdgeCoords()
    fea = MotorProblem(i_abc=i_abc, old_edge_coords=old_edge_coords)
    sim = Simulator(M(fea=fea))
    sim['edge_deltas'] = edge_deltas
    sim.run()
    plt.figure(1)
    fea.moveMesh(fea.uhat)
    plot(fea.mesh)
    plt.show()

#    print("CSDL: Running check_partials()...")
#    sim.check_partials()
