from fenics import *
import dolfin as do
import matplotlib.pyplot as plt
from msh2xdmf import import_mesh
import numpy as np
import time

set_log_level(1)

# Stator winding domains have been defined; these can have specific J values
# relating to current density in Maxwell's equations; 3-phase means phase shift of 2pi/3 rad
# between each current
# Boundary conditions set on x-axis and edge of domain
# individual domains have permeability depending on area (air, etc)
t_start = time.time()
tol = 1e-15
s = 9
p = 8

theta = 0
Imax = 282.8 # (A) Max current in wires defined in paper
# Imax = 94. # (A) No Load Current
Nstp = 13 # Number of turns per phase
# Nstp = 1 

theta_t = pi/s
theta_sso = .5 * theta_t
Rsy = 103.
Rssi = 83.

theta_p = 2*pi/2/p # Angular sweep of each Rotor Slice
theta_m = .78 * theta_p # Angular sweep of magnet
theta_b = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
theta_g = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet side piece

Rtm = 79.
Rbm = 74.5

Hc = 838.e3 # Coercivity of Neodymium magnet in kA/m (ranges between 800 and 950)
vacuum_perm = 4.e-7 * np.pi

winding_area = theta_sso/2/2 * (Rsy**2 - Rssi**2) * 10**(-6)
magnet_area = (Rtm**2 - Rbm**2) * theta_m/2

# Setting geometrical conditions for boundaries
class Domain_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        # return on_boundary and x[1] > DOLFIN_EPS
        return on_boundary

# Setting geometry for periodic boundary condition
class Periodic_Boundary(do.SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and x[0] > DOLFIN_EPS and near(x[1],0.)
        # Taking right side of semicircle boundary as original "target domain"

    def map(self,x,y): # Mapping to new boundary, w/ same conditions for negative part of x-axis  
        y[0] = -x[0]
        y[1] = x[1]

domain_bound = Domain_Boundary() # Far-Field Boundary Condition (approaches zero)
semicircle_bound = Periodic_Boundary() # Periodic Boundary Conditions
# -------------------------------------------------------------------------------------------
mesh, boundaries_mf, subdomains_mf, association_table = import_mesh(
    prefix="motor_mesh_new",
    dim=2,
    subdomains=True
)
# boundaries_mf: Mesh Function object for boundaries
# subdomains_mf: Mesh Function object for subdomains

V = FunctionSpace(mesh,'P',1, constrained_domain = Periodic_Boundary())

dx = Measure('dx', domain = mesh, subdomain_data = subdomains_mf)
dS = Measure('dS', domain = mesh, subdomain_data = boundaries_mf)


'''
----- Domain list for dx: -----
Rotor Core: 1
Stator Core: 2
Right Winding: 3 to 19 in intervals of 2
Left Winding: 4 to 20 in intervals of 2
Magnets: 21 to 28
Right Air Slots: 29 to 43 in intervals of 2
Left Air Slots: 30 to 44 in intervals of 2
Rotor Air Gap: 45
Middle Air Gap: 46 
Outer Domain: 47

----- Boundary list for ds: -----
Right Magnet Edge: 1 to 15 in intervals of 2
Left Magnet Edge: 2 to 16 in intervals of 2
'''

# Boundary Conditions
A_FF = Constant(0) # Far Field Magnetic Vector Potential (decays to zero)
bc0 = DirichletBC(V,A_FF,Domain_Boundary())

# Current Densities and Winding Indices
J1_value = Imax * Nstp * cos(theta) # Current Density 1 value
J2_value = Imax * Nstp * cos(theta + 2*pi/3)# Current Density 2 value
J3_value = Imax * Nstp * cos(theta - 2*pi/3) # Current Density 3 value

J1 = Constant(J1_value) # Current Density 1
J2 = Constant(J2_value) # Current Density 2
J3 = Constant(J3_value) # Current Density 3

WA_N = [3, 6, 7] # Winding Domain indices for current A North
WA_S = [4, 5, 8] # Winding Domain indices for current A South
WB_N = [9, 12, 13] # Winding Domain indices for current B North
WB_S = [10, 11, 14] # Winding Domain indices for current B South
WC_N = [15, 18, 19] # Winding Domain indices for current C North
WC_S = [16, 17, 20] # Winding Domain indices for current C South

# Define Magnetic Permeability
from piecewise_permeability import *

def RelativePermeability(subdomain, A_z):
    if subdomain == 1 or subdomain == 2:
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
        mu = 1000

    elif subdomain >= 3 and subdomain <= 20:
        mu = 1.0

    elif subdomain >= 21 and subdomain <= 28:
        mu = 1.05

    elif subdomain >= 29:
        mu = 1.0

    return mu

# Define Variational Problem
A_z_ = TrialFunction(V)
A_z = Function(V)
v = TestFunction(V)

# Current Densities for each winding phase
JA_N = sum(J1*v*dx(WA_N[i]) for i in range(3))
JA_S = -sum(J1*v*dx(WA_S[i]) for i in range(3))
JB_N = sum(J2*v*dx(WB_N[i]) for i in range(3))
JB_S = -sum(J2*v*dx(WB_S[i]) for i in range(3))
JC_N = sum(J3*v*dx(WC_N[i]) for i in range(3))
JC_S = -sum(J3*v*dx(WC_S[i]) for i in range(3))

# Current Densities for magnets
JM = 0.
num_magnets = 8
for i in range(num_magnets):
    angle = np.pi/num_magnets/2. + np.pi/num_magnets * i
    Hx = do.Constant((-1)**i * Hc * np.cos(angle))
    Hy = do.Constant((-1)**i * Hc * np.sin(angle))

    H = do.as_vector([Hx, Hy])
    curl_v = do.as_vector([v.dx(1), -v.dx(0)])

    JM += do.inner(H, curl_v) * dx(i+21)

J_W = JA_N + JA_S + JB_N + JB_S + JC_N + JC_S

L = JM + J_W

a = 0.
for i in range(46):
    a += 1. / (vacuum_perm) * 1/RelativePermeability(i + 1, A_z) * do.dot(do.grad(A_z), do.grad(v)) * dx(i + 1)

F = a - L
# A_z = do.Function(V)
# F = do.action(F, A_z)
Jac = derivative(F, A_z, A_z_)

problem = do.NonlinearVariationalProblem(F, A_z, bc0, Jac)
solver = do.NonlinearVariationalSolver(problem)
solver.parameters['nonlinear_solver']='snes' 
solver.parameters['snes_solver']['line_search'] = 'bt' 
solver.solve()

# solve(a - L == 0., A_z, bc0) # nonlinear solver

# Computing magnetic field (B = curl(A))
W = VectorFunctionSpace(mesh,'P',1)
B = project(as_vector((A_z.dx(1), -A_z.dx(0))),W)

t_end = time.time() - t_start

print('computation time:', t_end)

'''MESH PLOT'''
plt.figure(1)
# fig = plt.figure(1,figsize = (9,8))
# ax = fig.gca(projection="3d")
# ax.view_init(elev=90., azim=-90.)
# plt.axis([-4.5,4.5,0,4.25])
# do.plot(subdomains_mf)
# do.plot(A_z)
# do.plot(mesh, linewidth = .1)
do.plot(B, linewidth = 20)

plt.figure(2)
do.plot(A_z)

vtkfile_A_z = File('solutions/Magnetic_Vector_Potential.pvd')
vtkfile_B = File('solutions/Magnetic_Flux_Density.pvd')
vtkfile_A_z << A_z
vtkfile_B << B

plt.show()

# -------------------------------------------------------------------------------------------
'''TO-DO'''
# Set proper permeability values (and figure out the thing w vacuum)
# Assign the whole boundary at the x-axis 

# Figure out what the debug code means in output:
# DEBUG: [at /Users/lucascotzniovsky/opt/anaconda3/include/dolfin/mesh/MeshFunction.h:485 in operator=()]
# DEBUG: Mesh value collection does not contain all values for all entities

# for dx, we can call dx(n), where n is the index of the domain number in the Mesh Function
# we need this for the stator windings to assign the values of J1,J2,J3 to each winding
# This has to be done manually, but we can arrange a loop for currents in and out of the 