from dolfin import *
from csdl import CustomImplicitOperation, Model
import csdl
import numpy as np
from csdl_om import Simulator
from motor_fea import *


class MeshMotionModel(Model):

    def initialize(self):
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
        self.parameters.declare('fea')

    def define(self):
        self.fea = self.parameters['fea']
        self.input_size = len(self.fea.edge_indices)
        self.output_size = self.fea.total_dofs_uhat
        self.add_input('edge_deltas',
                        shape=(self.input_size,),
                        val=0.1*np.ones(self.input_size).reshape(self.input_size,)
                        )
        self.add_output('uhat',
                        shape=(self.output_size,),
                        val=0.0
                        )
        self.declare_derivatives('uhat', 'uhat')
        self.declare_derivatives('uhat', 'edge_deltas')

    def updateMeshMotionBC(self, inputs):
        self.fea.edge_deltas = inputs['edge_deltas']
        for i in range(self.input_size):
            self.fea.uhat_old.vector()[self.fea.edge_indices[i]] = inputs['edge_deltas'][i]
        self.bcs = self.fea.setBCMeshMotion()

    def evaluate_residuals(self, inputs, outputs, residuals):
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])
        resM = assemble(self.fea.resM())
        uhat_0 = Function(self.fea.VHAT)
        uhat_0.assign(self.fea.uhat)

        self.bcs.apply(resM, uhat_0.vector())
        residuals['uhat'] = resM.get_local()

    def solve_residual_equations(self, inputs, outputs):
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])

        self.fea.solveMeshMotion()

        outputs['uhat'] = self.fea.uhat.vector().get_local()
        update(self.fea.uhat, outputs['uhat'])


    def compute_derivatives(self, inputs, outputs, derivatives):
        self.updateMeshMotionBC(inputs)
        update(self.fea.uhat, outputs['uhat'])

        self.dRdu = assemble(self.fea.dRm_duhat)
        self.bcs.apply(self.dRdu)
        self.dRdf = self.fea.getBCDerivatives()
        self.A = self.dRdu


    def compute_jacvec_product(self, inputs, outputs,
                                d_inputs, d_outputs, d_residuals, mode):
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
    f = open('coarse_mesh_Ru/init_edge_coords.txt', 'r+')
    old_edge_coords = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    f = open('coarse_mesh_Ru/edge_coord_deltas.txt', 'r+')
    edge_deltas = np.fromstring(f.read(), dtype=float, sep=' ')
    f.close()

    print("number of nonzero displacements:", np.count_nonzero(edge_deltas))
    
    # One-time computation for the initial edge coordinates from
    # the code that creates the mesh file
    # old_edge_coords = getInitialEdgeCoords()
    
#    fea = MotorFEA(mesh_file="mesh_files/motor_mesh_1",
    fea = MotorFEA(mesh_file="coarse_mesh_Ru/motor_mesh_coarse_1",
                            old_edge_coords=old_edge_coords)
    sim = Simulator(MeshMotionModel(fea=fea))
    sim['edge_deltas'] = 0.5*edge_deltas
    sim.run()
    plt.figure(1)
    fea.moveMesh(fea.uhat)
    plot(fea.mesh)
    plt.show()
    
#Ru: check_partials may take a long running time; 
#tested to be correct

#    print("CSDL: Running check_partials()...")
#    sim.check_partials()
