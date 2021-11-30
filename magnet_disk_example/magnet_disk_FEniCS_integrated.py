import dolfin as do
import numpy as np
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh

do.set_log_level(1)

# PROBLEM SPECIFIC PARAMETERS
Hc = 838.e3  # 838 kA/m
num_magnets = 4
vacuum_perm = 4e-7 * np.pi

# MESH IMPORT
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh(
    prefix="magnet_disk_mesh",
    dim=2,
    subdomains=True
)

'''
SUBDOMAIN DATA:
DISK: 1
MAGNETS 1 - 4: 2 - 5

BOUNDARY DATA:
MAGNET IN-CURRENT LINES: ODD # UP TO 8
MAGNET OUT-CURRENT LINES: EVEN # UP TO 8
'''

# SETTING CLASSES FOR BOUNDARY CONDITIONS
class OuterBoundary(do.SubDomain):
    def inside(self, x, on_boundary):
        # return on_boundary and np.sqrt(x[0]**2 + x[1]**2) < Ri + do.DOLFIN_EPS
        return on_boundary
        # return on_boundary and np.sqrt(x[0]**2 + x[1]**2) > Ro - do.DOLFIN_EPS

outer_bound         = OuterBoundary()

V       = do.FunctionSpace(mesh, 'P', 1)
dx      = do.Measure('dx', domain=mesh, subdomain_data=subdomains_mf)
dS      = do.Measure('dS', domain=mesh, subdomain_data=boundaries_mf)               
A_outer = do.Constant(0)
bc_o    = do.DirichletBC(V, A_outer, outer_bound)

from piecewise_permeability import *

# START NEW PERMEABILITY
def RelativePermeability(subdomain, u, uhat):
    gradu = gradx(u,uhat)
    if subdomain == 1: # Electrical/Silicon/Laminated Steel
        B = do.as_vector((gradu[1], -gradu[0]))
        norm_B = do.sqrt(do.dot(B, B) + do.DOLFIN_EPS)
        
        mu = do.conditional(
            do.lt(norm_B, 1.004),
            linearPortion(norm_B),
            do.conditional(
                do.lt(norm_B, 1.433),
                cubicPortion(norm_B),
                (expA * do.exp(expB*norm_B + expC) + 1)
            )
        )
        # mu = 4000. # constant case

    else: # Magnets
        mu = 1.05

    return mu

# END NEW PERMEABILITY

VHAT = do.VectorFunctionSpace(mesh, 'P', 1)
uhat = do.Function(VHAT)

A_z     = do.Function(V)
v       = do.TestFunction(V)


I = do.Identity(2)
def gradx(f,u):
    return do.dot(do.grad(f), do.inv(I + do.grad(u)))

def J(uhat):
    return do.det(I + do.grad(uhat))

def Jm(v,uhat,num_magnets,Hc):
    Jm = 0.
    gradv = gradx(v,uhat)
    for i in range(num_magnets):
        angle = np.pi / 4 + i * np.pi / 2
        Hx = do.Constant((-1)**i * Hc * np.cos(angle))
        Hy = do.Constant((-1)**i * Hc * np.sin(angle))

        H = do.as_vector([Hx, Hy])
        
        curl_v = do.as_vector([gradv[1],-gradv[0]])
        Jm += do.inner(H,curl_v)*dx(i+2)
    return Jm
    
def pdeRes(u,v,uhat,num_magnets,Hc):
    res = 0
    gradu = gradx(u,uhat)
    gradv = gradx(v,uhat)
    for i in range(5):
        res += 1./vacuum_perm*(1/RelativePermeability(i + 1, u, uhat))\
                *do.dot(gradu,gradv)*dx(i + 1)
    res -= Jm(v,uhat,num_magnets,Hc)
    return res

######### Subroutine to move the mesh #################################
from magnet_disk_mesh import generateMeshMovement

#_, old_edge_coords, edge_deltas, edge_indices = generateMeshMovement()
# We will then manually find the `edge_indices` in FEniCS by KDTree
_, old_edge_coords, edge_deltas, _ = generateMeshMovement()
uhat_0 = do.Function(VHAT)

# TODO: need to make sure the numbering of these two are the same
V0 = do.FunctionSpace(mesh, 'CG', 1)
coordinates = V0.tabulate_dof_coordinates()
from scipy.spatial import KDTree
def findNodeIndices(node_coordinates, coordinates):
    tree = KDTree(coordinates)
    dist, node_indices = tree.query(node_coordinates)
    return node_indices
    
node_indices = findNodeIndices(np.reshape(old_edge_coords, (-1,2)), coordinates)
edge_indices = np.empty(2*len(node_indices))
for i in range(len(node_indices)):
    edge_indices[2*i] = 2*node_indices[i]
    edge_indices[2*i+1] = 2*node_indices[i]+1
    
# Verifying the accuracy of the indices
#print(np.linalg.norm(coordinates[node_indices]-np.reshape(old_edge_coords, (-1,2)), axis=1))

# Assigning the displacements to corresponding dofs of `uhat_0`
uhat_0.vector()[edge_indices] = edge_deltas

# A facet MeshFunction of dim=mesh.topoloty().dim()-1 is needed to set the BC
bc_moved = do.DirichletBC(VHAT, uhat_0, boundaries_mf, 1000)

vhat = do.TrialFunction(VHAT)
a = do.inner(do.grad(uhat), do.grad(vhat))*do.dx
L = 0
do.solve(a==L, uhat, bcs=bc_moved)
plt.figure(1)
do.plot(uhat)
do.ALE.move(mesh, uhat)
plt.figure(2)
do.plot(mesh)
vtkfile_mesh = do.File('magnet_disk_solutions/mesh_movement.pvd')
vtkfile_mesh << uhat
plt.show()
#######################################################################

#pdeRes = pdeRes(A_z,v,uhat,num_magnets,Hc)
#J = do.derivative(pdeRes,A_z)
#do.solve(pdeRes==0,A_z,J=J,bcs=bc_o)
##do.ALE.move(mesh, uhat)
#W       = do.VectorFunctionSpace(mesh,'DG',0)
#B       = do.project(do.as_vector((A_z.dx(1),-A_z.dx(0))),W)

#plt.figure(1)
#do.plot(B,linewidth=40)

#plt.figure(2)
#do.plot(A_z)

#plt.figure(3)
#do.plot(mesh)

#vtkfile_A_z = do.File('magnet_disk_solutions/Magnetic_Vector_Potential.pvd')
#vtkfile_B = do.File('magnet_disk_solutions/Magnetic_Flux_Density.pvd')
#vtkfile_A_z << A_z
#vtkfile_B << B

#plt.show()
