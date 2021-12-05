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
    prefix="magnet_disk_mesh_1",
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
def RelativePermeability(subdomain, A_z):
    if subdomain == 1: # Electrical/Silicon/Laminated Steel
        B = do.as_vector((A_z.dx(1), -A_z.dx(0)))
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

A_z     = do.Function(V)
v       = do.TestFunction(V)

# A_z as Function
# setting up nonlinear  problem using residual (r == a-L == 0 vs. a == L)

# SETTING UP EQUIVALENT CURRENTS IN MAGNETS
Jm = 0.
for i in range(num_magnets):
    angle = np.pi / 4 + i * np.pi / 2
    Hx = do.Constant((-1)**i * Hc * np.cos(angle))
    Hy = do.Constant((-1)**i * Hc * np.sin(angle))

    H = do.as_vector([Hx, Hy])
    curl_v = do.as_vector([v.dx(1), -v.dx(0)])  
    Jm += do.inner(H,curl_v)*dx(i+2)

# a       = (1/mu)*do.dot(do.grad(A_z),do.grad(v))*dx
a = 0
for i in range(5):
    a += 1./vacuum_perm*(1/RelativePermeability(i + 1, A_z))*do.dot(do.grad(A_z),do.grad(v))*dx(i + 1)

L       = Jm

# do.solve(a==L,A_z,bc_o) # linear solver
do.solve(a-L == 0.,A_z,bc_o) # nonlinear solver

W       = do.VectorFunctionSpace(mesh,'DG',0)
B       = do.project(do.as_vector((A_z.dx(1),-A_z.dx(0))),W)

plt.figure(1)
do.plot(B, linewidth=40)

plt.figure(2)
do.plot(A_z)

vtkfile_A_z = do.File('magnet_disk_solutions/Magnetic_Vector_Potential.pvd')
vtkfile_B = do.File('magnet_disk_solutions/Magnetic_Flux_Density.pvd')
vtkfile_A_z << A_z
vtkfile_B << B

plt.show()