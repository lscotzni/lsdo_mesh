import lsdo_mesh.lsdo_mesh as lm
import numpy as np

'''
Mesh Generation of 2D Radial Flux PMSM; only half is modeled (180 degrees) due to periodicity
'''

p       = 8 # pole pairs
s       = 9 # number of stator slots for windings per 180 degrees
m       = 3 # number of phases for stator winding current

ks      = 1 # target mesh size

''' -------------------- Key Geometrical Parameters of Motor -------------------- '''

# # Key Rotor Angles for Geometry
theta_p = 2*pi/2/p # Angular sweep of each Rotor Slice
theta_m = .78 * theta_p # Angular sweep of magnet
theta_b = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
theta_g = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet side piece

# # Key Stator Angles for Geometry
theta_t = pi/s # Angular sweep of each Stator Tooth
theta_sso = .5 * theta_t # Angular sweep of total area of windings
theta_ssi = .3 * theta_sso # Angular sweep of tooth tip separation

# Need to define each of these
Rr = 80.
Rtm = 79.
Rtb = 77.5
Rbb = 75.
Rbm = 74.5
Rin = 60. # Inner Radius of Rotor

Rout = 115.
Rsy = 103.
Rssi = 83.
Rs = 81.

# D = 200. # Domain Radius
D = 175. # Domain Radius
Ras = (Rr+Rs)/2 # Radius splitting air-gap mesh of rotor and stator
RR = (Rin+Rr)/2 # Midpoint to cut air-gap mesh in Rotor
RS = (Rsy+Rout)/2 # # Midpoint to cut air-gap mesh in Stator

