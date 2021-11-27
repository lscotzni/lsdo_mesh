from lsdo_mesh import geometry_classes as lm
import numpy as np
import matplotlib.pyplot as plt

'''
Mesh Generation of 2D Radial Flux PMSM; only half is modeled (180 degrees) due to periodicity
'''
''' -------------------- Motor Attributes -------------------- '''
p       = 8 # pole pairs
s       = 9 # number of stator slots for windings per 180 degrees
m       = 3 # number of phases for stator winding current

ks      = 10 # target mesh size

''' -------------------- Key Geometrical Parameters of Motor -------------------- '''
# # Key Rotor Angles for Geometry
theta_p     = 2*np.pi/2/p # Angular sweep of each Rotor Slice
theta_m     = .78 * theta_p # Angular sweep of magnet
theta_b     = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
theta_g     = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet side piece

# # Key Stator Angles for Geometry
theta_t     = np.pi/s # Angular sweep of each Stator Tooth
theta_sso   = .5 * theta_t # Angular sweep of total area of windings
theta_ssi   = .3 * theta_sso # Angular sweep of tooth tip separation

# Need to define each of these
Rr          = 80.
Rtm         = 79.
Rtb         = 77.5
Rbb         = 75.
Rbm         = 74.5
Rin         = 60. # Inner Radius of Rotor

Rout        = 115.
Rsy         = 103.
Rssi        = 83.
Rs          = 81.

# D = 200. # Domain Radius
D           = 175. # Domain Radius
Ras         = (Rr+Rs)/2 # Radius splitting air-gap mesh of rotor and stator
RR          = (Rin+Rr)/2 # Midpoint to cut air-gap mesh in Rotor
RS          = (Rsy+Rout)/2 # # Midpoint to cut air-gap mesh in Stator

m           = lm.Mesh(name='motor_mesh_new', popup=False)

# NOTE: need to fix implementation; the process does not like lists
# as inputs, but the option to do so makes things easier (especially 
# this motor model process)
#   - NOTE: This has been done for lm.Surface, still need to apply for others

# ---------------------------- POINTS ----------------------------
origin      = lm.Point(0., 0. ,ms=ks)

# -- Inner Rotor
rotor_inner_surface_p   = [lm.Point(Rin, np.pi * i / 2, ms=ks, mode='polar') for i in range(3)]

# -- Outer Rotor
rotor_outer_surface_p = [
    lm.Point(Rr, np.pi / p * i, ms=ks, mode='polar') for i in range(p+1)
] # not sure if it's better to define per slice or just as 2 circle arcs

# -- Magnets and Air Slots
magnet_air_slot_1_p     = [
    lm.Point(Rbm, theta_b / 2 + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_2_p     = [
    lm.Point(Rtm, theta_b / 2 + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_3_p     = [
    lm.Point(Rbm, theta_b / 2 + theta_g + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_4_p     = [
    lm.Point(Rtm, theta_b / 2 + theta_g + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_5_p     = [
    lm.Point(Rbm, theta_b / 2 + theta_g + theta_m + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_6_p     = [
    lm.Point(Rtm, theta_b / 2 + theta_g + theta_m + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_7_p     = [
    lm.Point(Rbm, theta_p - theta_b / 2 + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

magnet_air_slot_8_p    = [
    lm.Point(Rtm, theta_p - theta_b / 2 + np.pi / p * i, ms=ks, mode='polar') for i in range(p)
]

# -- Stator Outer Surface
# stator_outer_surface_p  = [lm.Point(Rout, 0, ms=ks)]
# [stator_outer_surface_p.extend(lm.rotate([stator_outer_surface_p[0]], angle=[0., 0., (i + 1) * np.pi / 2])) for i in range(2)]

stator_outer_surface_p = [
    lm.Point(Rout, np.pi * i / 2, ms=ks, mode='polar') for i in range(3)
]

# -- Stator Inner Surface
stator_inner_surface_1_p = [
    lm.Point(Rsy, np.pi / s * i, ms=ks, mode='polar') for i in range(s + 1)
]

stator_inner_surface_2_p = [
    lm.Point(Rsy, theta_sso / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_3_p = [
    lm.Point(Rssi, theta_sso / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_4_p = [
    lm.Point(Rssi, theta_ssi / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_5_p = [
    lm.Point(Rs, theta_ssi / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_6_p = [
    lm.Point(Rs, theta_t - theta_ssi / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_7_p = [
    lm.Point(Rssi, theta_t - theta_ssi / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_8_p = [
    lm.Point(Rssi, theta_t - theta_sso / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

stator_inner_surface_9_p = [
    lm.Point(Rsy, theta_t - theta_sso / 2 + np.pi / s * i, ms=ks, mode='polar') for i in range(s)
]

# -- Domain Boundary
domain_boundary_p       = [lm.Point(D, 0, ms=ks)]
for i in range(2):
    domain_boundary_p.extend(lm.rotate([domain_boundary_p[0]], angle = [0., 0., (i + 1) * np.pi / 2]))

# -- Stator Winding Mid-Section
stator_winding_mid_p    = [lm.Point(Rssi, 0., ms=ks)]
[stator_winding_mid_p.extend(lm.rotate([stator_winding_mid_p[0]], angle = [0., 0., (i + 1) * np.pi / s])) for i in range(s - 1 + 1)]

# -- Stator Core Circle Sector Points
stator_core_circle_sector_p = [
    lm.Point((Rr) / 2., 0.,  ms=ks, mode='polar'),
    lm.Point(Rout, 0., ms=ks, mode='polar'),
    lm.Point(Rout, np.pi / 2., ms=ks, mode='polar'),
    lm.Point(Rout, np.pi, ms=ks, mode='polar'),
    lm.Point((Rr) / 2., np.pi / 2.,  ms=ks, mode='polar'),
    lm.Point((Rr) / 2., np.pi,  ms=ks, mode='polar')
]

# ---------------------------- CURVES ----------------------------
# -- Rotor Inner Surface
rotor_inner_surface_c   = [
    lm.Curve(rotor_inner_surface_p[0] ,rotor_inner_surface_p[1] ,origin ,curve_type = 'arc'),
    lm.Curve(rotor_inner_surface_p[1] ,rotor_inner_surface_p[2] ,origin ,curve_type = 'arc'),
    lm.Curve(rotor_inner_surface_p[2] ,rotor_inner_surface_p[0]),
]

# -- Rotor Outer Surface
rotor_outer_surface_c   = []
for i in range(p):
    rotor_outer_surface_c.append(
        lm.Curve(rotor_outer_surface_p[-1 - i], rotor_outer_surface_p[-2 - i], origin, curve_type='arc')
    )

rotor_outer_surface_c.append(lm.Curve(rotor_outer_surface_p[-1], rotor_outer_surface_p[0]))

# -- Magnet & Air Slot Curves
right_air_slot_c, left_air_slot_c = [], []
for i in range(p):
    right_air_slot_c.extend([[
        lm.Curve(magnet_air_slot_1_p[i], magnet_air_slot_2_p[i]),
        lm.Curve(magnet_air_slot_2_p[i], magnet_air_slot_4_p[i], origin, curve_type='arc'),
        lm.Curve(magnet_air_slot_4_p[i], magnet_air_slot_3_p[i]),
        lm.Curve(magnet_air_slot_3_p[i], magnet_air_slot_1_p[i], origin, curve_type='arc'),
    ]])
    left_air_slot_c.extend([[
        lm.Curve(magnet_air_slot_5_p[i], magnet_air_slot_6_p[i]),
        lm.Curve(magnet_air_slot_6_p[i], magnet_air_slot_8_p[i], origin, curve_type='arc'),
        lm.Curve(magnet_air_slot_8_p[i], magnet_air_slot_7_p[i]),
        lm.Curve(magnet_air_slot_7_p[i], magnet_air_slot_5_p[i], origin, curve_type='arc'),
    ]])

magnet_c = []
for i in range(p):
    magnet_c.extend([[
        lm.Curve(magnet_air_slot_3_p[i], magnet_air_slot_4_p[i], physical_group=(2*i+1, 'Magnet Side R '+str(i+1))),
        lm.Curve(magnet_air_slot_4_p[i], magnet_air_slot_6_p[i], origin, curve_type='arc'),
        lm.Curve(magnet_air_slot_6_p[i], magnet_air_slot_5_p[i], physical_group=(2*i+2, 'Magnet Side L '+str(i+1))),
        lm.Curve(magnet_air_slot_5_p[i], magnet_air_slot_3_p[i], origin, curve_type='arc'),
    ]])

right_air_slot_c = list(
    np.array(right_air_slot_c).reshape((p * len(right_air_slot_c[0])))
)
left_air_slot_c = list(
    np.array(left_air_slot_c).reshape((p * len(left_air_slot_c[0])))
)
magnet_c = list(
    np.array(magnet_c).reshape((p * len(magnet_c[0])))
)

# -- Outer Stator Surface Curves
outer_stator_surface_c  = [
    lm.Curve(stator_outer_surface_p[0], stator_outer_surface_p[1], origin, curve_type='arc'),
    lm.Curve(stator_outer_surface_p[1], stator_outer_surface_p[2], origin, curve_type='arc')
]

# -- Stator Surface Entites (Stator Core + Windings)
stator_core_boundary_c  = []
stator_core_boundary_c.extend([[
    lm.Curve(stator_inner_surface_1_p[i], stator_inner_surface_2_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_2_p[i], stator_inner_surface_3_p[i]),
    lm.Curve(stator_inner_surface_3_p[i], stator_inner_surface_4_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_4_p[i], stator_inner_surface_5_p[i]),
    lm.Curve(stator_inner_surface_5_p[i], stator_inner_surface_6_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_6_p[i], stator_inner_surface_7_p[i]),
    lm.Curve(stator_inner_surface_7_p[i], stator_inner_surface_8_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_8_p[i], stator_inner_surface_9_p[i]),
    lm.Curve(stator_inner_surface_9_p[i], stator_inner_surface_1_p[i + 1], origin, curve_type='arc'),
] for i in range(s)])

stator_core_boundary_c  = list(
    np.array(stator_core_boundary_c).reshape(s*len(stator_core_boundary_c[0]))
)
stator_core_boundary_c.extend([
    lm.Curve(stator_inner_surface_1_p[-1], stator_outer_surface_p[-1]),
    lm.Curve(stator_outer_surface_p[1], stator_outer_surface_p[2], origin, curve_type='arc'),
    lm.Curve(stator_outer_surface_p[0], stator_outer_surface_p[1], origin, curve_type='arc'),
    lm.Curve(stator_outer_surface_p[0], stator_inner_surface_1_p[0])
])

# -- Stator Winding Curves
stator_windings_right_c = []
stator_windings_right_c.extend([[
    lm.Curve(stator_inner_surface_1_p[i], stator_inner_surface_2_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_2_p[i], stator_inner_surface_3_p[i]),
    lm.Curve(stator_inner_surface_3_p[i], stator_inner_surface_4_p[i], origin,  curve_type='arc'),
    lm.Curve(stator_inner_surface_4_p[i], stator_winding_mid_p[i], origin,  curve_type='arc'),
    lm.Curve(stator_winding_mid_p[i], stator_inner_surface_1_p[i])
] for i in range(s)])

stator_windings_left_c = []
stator_windings_left_c.extend([[
    lm.Curve(stator_inner_surface_9_p[i], stator_inner_surface_1_p[i+1], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_1_p[i+1], stator_winding_mid_p[i+1]),
    lm.Curve(stator_winding_mid_p[i+1], stator_inner_surface_7_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_7_p[i], stator_inner_surface_8_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_8_p[i], stator_inner_surface_9_p[i]),
] for i in range(s)])

stator_windings_right_c = list(
    np.array(stator_windings_right_c).reshape((s * len(stator_windings_right_c[0])))
)

stator_windings_left_c = list(
    np.array(stator_windings_left_c).reshape((s * len(stator_windings_left_c[0])))
)

# -- Air Gap Curves
air_gap_curves = []
air_gap_curves.extend([[
    lm.Curve(stator_winding_mid_p[i], stator_inner_surface_4_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_4_p[i], stator_inner_surface_5_p[i], origin),
    lm.Curve(stator_inner_surface_5_p[i], stator_inner_surface_6_p[i], origin, curve_type='arc'),
    lm.Curve(stator_inner_surface_6_p[i], stator_inner_surface_7_p[i], origin),
    lm.Curve(stator_inner_surface_7_p[i], stator_winding_mid_p[i + 1], origin, curve_type='arc'),
] for i in range(s)])

air_gap_curves = list(
    np.array(air_gap_curves).reshape((s * len(air_gap_curves[0])))
)

# air_gap_curves.extend(stator_entities_boundary_c[:45])
air_gap_curves.append(lm.Curve(stator_winding_mid_p[-1],rotor_outer_surface_p[-1]))
air_gap_curves.extend([rotor_outer_surface_c[i] for i in range(8)])
air_gap_curves.append(lm.Curve(rotor_outer_surface_p[0], stator_winding_mid_p[0]))

# air_gap_curves.extend(rotor_outer_surface_c[:-1])
# air_gap_curves.append(lm.Curve(rotor_outer_surface_p[-1], stator_winding_mid_p[-1]))
# air_gap_curves.extend([stator_entities_boundary_c[-(i+5)] for i in range(45)])
# air_gap_curves.append(lm.Curve(stator_winding_mid_p[0], rotor_outer_surface_p[0]))

# -- Domain Boundary Curves
domain_boundary_c       = [
    lm.Curve(domain_boundary_p[0], domain_boundary_p[1], origin, curve_type='arc'),
    lm.Curve(domain_boundary_p[1], domain_boundary_p[2], origin, curve_type='arc'),
    lm.Curve(domain_boundary_p[2], stator_outer_surface_p[2]),
    # outer_stator_surface_c[1],
    # outer_stator_surface_c[0],
    lm.Curve(stator_outer_surface_p[2], stator_outer_surface_p[1], origin, curve_type='arc'),
    lm.Curve(stator_outer_surface_p[0], stator_outer_surface_p[1], origin, curve_type='arc'), # this line is causing problems
    lm.Curve(stator_outer_surface_p[0], domain_boundary_p[0])
]

# NOTE: line 236 when the indices are swapped is interpreted as the same curve right above it, which is incorrect
#   - might be an issue with the duplicate removal but unsure

# ---------------------------- SURFACES ----------------------------
rotor_inner_surface_s   = lm.Surface(
    [rotor_inner_surface_c[0], rotor_inner_surface_c[1], rotor_inner_surface_c[2]],
    physical_group=(45, "Rotor Air Gap")
)
# m.add_entity(rotor_inner_surface_s)

# Magnets and corresponding air slots
magnet_s, left_air_slot_s, right_air_slot_s = [], [], []
for i in range(p):
    magnet_s.append(lm.Surface(magnet_c[4 * i:4 * i + 4], physical_group=(20 + i + 1, 'Magnet '+str(i+1))))
    left_air_slot_s.append(lm.Surface(left_air_slot_c[4 * i:4 * i + 4], physical_group=(20 + p + 2*i + 1, 'Right Air Slot '+str(i+1))))
    right_air_slot_s.append(lm.Surface(right_air_slot_c[4 * i:4 * i + 4], physical_group=(20 + p + 2*i + 2, 'Left Air Slot '+str(i+1))))

rotor_core_surface = []
rotor_core_surface.extend(rotor_outer_surface_c)

rotor_core_surface      = lm.Surface(rotor_core_surface)

rotor_core_subtract = []
rotor_core_subtract.append(rotor_inner_surface_s)
rotor_core_subtract.extend(magnet_s)
rotor_core_subtract.extend(left_air_slot_s)
rotor_core_subtract.extend(right_air_slot_s)
rotor_core              = lm.BooleanSurface(
    [rotor_core_surface], rotor_core_subtract, removeTool=False, operation='subtract', physical_group=(1, 'Rotor Core')
)
# rotor_core = rotor_core_surface
m.add_entity(rotor_core)

stator_core_surface     = lm.Surface(stator_core_boundary_c, physical_group=(2, 'Stator Core'))
m.add_entity(stator_core_surface)

right_winding_surfaces, left_winding_surfaces = [], []
for i in range(s):
    right_winding_surfaces.append(lm.Surface(stator_windings_right_c[5  * i: 5 * i + 5], physical_group=(2+(2*i+1), 'Right Winding '+str(i+1))))
    left_winding_surfaces.append(lm.Surface(stator_windings_left_c[5  * i: 5 * i + 5], physical_group=(2+(2*i+2), 'Left Winding '+str(i+1))))

    m.add_entity(right_winding_surfaces[i])
    m.add_entity(left_winding_surfaces[i])

# domain_boundary_surface     = lm.Surface(domain_boundary_c, physical_group = (47, 'Outer Domain'))
# m.add_entity(domain_boundary_surface)

air_gap_surface = lm.Surface(air_gap_curves, physical_group = (46, 'Air Gap'))
m.add_entity(air_gap_surface)

# m.add_all_entities_to_physical_group(geometry_type='curves')



# FFD PORTION

mesh_points  = m.get_coordinates(coord_sys='cartesian')
# print(mesh_points)

plt.plot(
    [point[0] for point in mesh_points],
    [point[1] for point in mesh_points],
    'k*'
)
# plt.show()
inner_rotor_f = []
inner_rotor_f.extend([
    lm.Face(
        rotor_inner_surface_p[i], 
        input_type='polar'
    ) for i in range(len(rotor_inner_surface_p))
])
inner_rotor_f[0].add_shape_parameter(
    'Inner Rotor Radius', 
    'r', 
    'constant'
)

m.assemble(coordinate_system='polar')

# inner_rotor_f[0].return_coordinates()


# [inner_rotor_f[i].add_shape_parameter(
#     'Inner Rotor Radius', 
#     [-0.5, 0.5], 
#     'r', 
#     1
# ) for i in range(len(rotor_inner_surface_p))]
# [print(inner_rotor_f[i].return_coordinates()) for i in range(len(rotor_inner_surface_p))]


# NOTE-S:
'''
Incomplete tasks:
- Curves for stator windings, magnets (& air gaps) and rotor-stator air gap
- Surfaces for stator windings, magnets (& air gaps) and rotor-stator air gap
- Boolean operations for stator windings, magnets (& air gaps) and rotor-stator air gap

ERRORS:
- if we try to apply boolean subtraction with a coincident boundary after a geometry 
    has been removed, the curve technically doesn't exist or isn't read so we cannot
    use it anymore; thus, stator, windings and air gaps need to be made manually
- we can make the windings and air gap, and then use these as tools to remove from a circle
    sector to create the stator core
'''


# os.system('python3 msh2xdmf.py -d 2 motor_mesh_new.msh')