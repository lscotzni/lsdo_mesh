import lsdo_mesh as lm
from math import pi
import numpy as np

'''
MOTOR MESH
'''

''' ========================== GEOMETRICAL PARAMETERS ========================== '''
p = 8
s = 9
m = 3

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

''' ========================== SET UP POINTS & MESH OBJECT ========================== '''
m   = lm.MotorMesh("Motor",pole_pairs = p, stator_teeth = s) # CALL TO INITALIZE GMSH

ms  = 1 # Mesh size

# Rotor End Point (prevent duplicate); we want to remove this later & automate within lsdo_mesh
outer_rotor_p_e     = lm.Point(-Rtm,0.0,ms = ms)

# Outer Surface of Rotor Core
outer_rotor_p_1     = lm.Point(Rtm,0.0,ms = ms)
outer_rotor_p_2     = lm.Point(Rtm,theta_b/2,ms = ms, coordinates = 'polar')
outer_rotor_p_3     = lm.Point(Rr, theta_b/2+theta_g, ms = ms, coordinates = 'polar')
outer_rotor_p_4     = lm.Point(Rr, theta_b/2+theta_g+theta_m, ms = ms, coordinates = 'polar')
outer_rotor_p_5     = lm.Point(Rtm, theta_p-theta_b/2, ms = ms, coordinates = 'polar')

# Magnets and Air Slots
slot_right_bottom_r = lm.Point(Rbb,theta_b/2, ms = ms, coordinates = 'polar')
slot_right_top_r    = lm.Point(Rtb,theta_b/2, ms = ms, coordinates = 'polar')
magnet_bottom_r     = lm.Point(Rbm,theta_b/2 + theta_g, ms = ms, coordinates = 'polar')
slot_right_bottom_l = lm.Point(Rbb,theta_b/2 + theta_g, ms = ms, coordinates = 'polar')
magnet_top_r        = lm.Point(Rtm,theta_b/2 + theta_g, ms = ms, coordinates = 'polar')
magnet_bottom_l     = lm.Point(Rbm,theta_b/2 + theta_g + theta_m, ms = ms, coordinates = 'polar')
slot_left_bottom_r  = lm.Point(Rbb,theta_b/2 + theta_g + theta_m, ms = ms, coordinates = 'polar')
magnet_top_l        = lm.Point(Rtm,theta_b/2 + theta_g + theta_m, ms = ms, coordinates = 'polar')
slot_left_bottom_l  = lm.Point(Rbb,theta_p - theta_b/2, ms = ms, coordinates = 'polar')
slot_left_top_l     = lm.Point(Rtb,theta_p - theta_b/2, ms = ms, coordinates = 'polar')

# Stator End Points (to prevent duplicates); we want to remove this and automate it later
inner_stator_p_e    = lm.Point(-Rsy, 0.0, ms = ms)

# Inner Surface of Stator Core & Winding Definition
inner_stator_p_1    = lm.Point(Rsy,0.0,ms = ms, coordinates = 'polar')
inner_stator_p_2    = lm.Point(Rsy,theta_sso/2,ms = ms, coordinates = 'polar')
inner_stator_p_3    = lm.Point(Rssi,theta_sso/2,ms = ms, coordinates = 'polar')
inner_stator_p_4    = lm.Point(Rssi,theta_ssi/2,ms = ms, coordinates = 'polar')
inner_stator_p_5    = lm.Point(Rs,theta_ssi/2,ms = ms, coordinates = 'polar')
inner_stator_p_6    = lm.Point(Rs,theta_t - theta_ssi/2,ms = ms, coordinates = 'polar')
inner_stator_p_7    = lm.Point(Rssi,theta_t - theta_ssi/2,ms = ms, coordinates = 'polar')
inner_stator_p_8    = lm.Point(Rssi,theta_t - theta_sso/2,ms = ms, coordinates = 'polar')
inner_stator_p_9    = lm.Point(Rsy,theta_t - theta_sso/2,ms = ms, coordinates = 'polar')

# Winding Corner in Mid-Section
winding_corner_p    = lm.Point(Rssi,0,ms = ms, coordinates = 'polar')

# Stator Winding Points (our case would be duplicates, but in general they may not be)


# ------------------------------ OBJECT TESTING START ------------------------------
c1 = lm.Curve(inner_stator_p_1,inner_stator_p_2,curve_type = 'line')
c2 = lm.Curve(inner_stator_p_1,inner_stator_p_2,inner_stator_p_3,curve_type = 'arc')

# print(c1)
# print(vars(c1))

# print(c2)
# print(vars(c2))

asdf  = lm.Point(np.sqrt(2),np.sqrt(2)).rotate(pi/4)
asdf2 = asdf.rotate(pi/4)

# print(vars(asdf))
# print(vars(asdf2))

exit()
# ------------------------------ OBJECT TESTING END ------------------------------


m.define_domain_boundary(D,ms)
m.define_shaft_surface(Rin,ms)
m.define_outer_surface(Rout,ms)



m.generate_mesh(points = points, ms = ms)

# print(vars(m))
print(vars(m)['rotor_points'])
# print(vars(m)['rotor_points'][0])
# print(vars(m)['rotor_points'][0].__dict__)


'''
NOTES:
- we can remove the points defining the inner radius of rotor (where it meets the shaft)
  and those defining the outer radius of the stator AND the domain points
    - we know where these points will always be as long we have the radius, so there isn't
      a need to manually define these points
'''


# pt1 = lm.Point(...)
# pt2 = lm.Point(...)
# cpt = lm.Point(...)
# c = lm.Curve()
# c2  = lm.Curve()
# s  = lm.Surface(c1,c2)
# m.add_surface(s)
# m.assemble() 