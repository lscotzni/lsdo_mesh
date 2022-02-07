from lsdo_mesh import geometry_classes as lm
import numpy as np
import matplotlib.pyplot as plt
import os

'''
Mesh Generation of 2D Radial Flux PMSM
'''
# =========================================================
#   !!!!!!! NOTE: NEED TO WRITE AS A FUNCTION !!!!!!!
# =========================================================

def MotorMeshGenerator(rotation_angles, file_name, poles):

    ''' -------------------- Motor Attributes -------------------- '''
    p       = poles # poles per 360 degrees
    s       = 3 * p # stator slots per 360 degrees
    m       = 3 # number of phases for stator winding current

    coarse_test = False
    if not coarse_test:
        # TARGET MESH SIZES (DIFFER FOR POINTS)
        l1      = 1e-2
        l2      = 3e-3
        l3      = 1e-3
        l4      = 3e-4
        ks      = 1e-2 # target mesh size
    elif coarse_test:
        l1      = 1e-2
        l2      = 1e-2
        l3      = 1e-2
        l4      = 1e-2
        ks      = 1e-2 # target mesh size

    # --- COLORS (r, g, b) ---
    white   = [255, 255, 255]
    red     = [255, 0, 0]
    blue    = [0, 0, 255]
    black   = [0, 0, 0]
    pink    = [255, 0, 255]
    green   = [0, 255, 0]
    yellow  = [255, 255, 0]

    ''' -------------------- Key Geometrical Parameters of Motor -------------------- '''
    # # Key Rotor Angles for Geometry
    theta_p     = 2*np.pi/p # Angular sweep of each Rotor Slice
    theta_m     = .6 * theta_p # Angular sweep of magnet
    # theta_b     = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
    # theta_g     = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet air gap
    theta_g     = .1 * theta_m

    # # Key Stator Angles for Geometry
    theta_t     = 2.0 * np.pi / s # Angular sweep of each Stator Tooth
    theta_sso   = .5 * theta_t # Angular sweep of total area of windings
    theta_ssi   = .3 * theta_sso # Angular sweep of tooth tip separation

    # Need to define each of these
    Rr          = 80.e-3
    Rtm         = 79.e-3
    Rtb         = 77.5e-3
    Rbb         = 75.e-3
    Rbm         = 74.5e-3
    # Rin         = 50.e-3 # Inner Radius of Rotor
    Rin         = 35.e-3

    Rout        = 120.e-3
    Rsy         = 103.e-3
    Rssi        = 83.e-3
    Rs          = 81.e-3

    D           = 175.e-3 # Domain Radius
    Ras         = (Rr+Rs)/2 # Radius splitting air-gap mesh of rotor and stator
    RR          = (Rin+Rr)/2 # Midpoint to cut air-gap mesh in Rotor
    RS          = (Rsy+Rout)/2 # # Midpoint to cut air-gap mesh in Stator

    rho         = rotation_angles
    m           = lm.Mesh(name=file_name, popup=True, rotation_angles=rho)

    # NOTE: need to fix implementation; the process does not like lists
    # as inputs, but the option to do so makes things easier (especially 
    # this motor model process)
    #   - NOTE: This has been done for lm.Surface, still need to apply for others

    # ---------------------------- POINTS ----------------------------
    origin      = lm.Point(0., 0., ms=l1)

    # -- Inner Rotor
    rotor_inner_surface_p   = [lm.Point(Rin, np.pi * i / 2, ms=l1, mode='polar') for i in range(4)]

    # -- Outer Rotor
    rotor_outer_surface_p = [
        lm.Point(Rr, np.pi / 2 * i, ms=l4, mode='polar') for i in range(4)
    ]

    # -- Magnets and Air Slots
    magnet_air_slot_1_p     = [
        # lm.Point(Rbm, theta_b / 2 + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rbm, 2*np.pi/p*i - theta_m/2 - theta_g, ms=l3, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_2_p     = [
        # lm.Point(Rtm, theta_b / 2 + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rtm, 2*np.pi/p*i - theta_m/2 - theta_g, ms=l4, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_3_p     = [
        # lm.Point(Rbm, theta_b / 2 + theta_g + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rbm, 2*np.pi/p*i - theta_m/2, ms=l3, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_4_p     = [
        # lm.Point(Rtm, theta_b / 2 + theta_g + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rtm, 2*np.pi/p*i - theta_m/2, ms=l4, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_5_p     = [
        # lm.Point(Rbm, theta_b / 2 + theta_g + theta_m + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rbm, 2*np.pi/p*i + theta_m/2, ms=l3, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_6_p     = [
        # lm.Point(Rtm, theta_b / 2 + theta_g + theta_m + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rtm, 2*np.pi/p*i + theta_m/2, ms=l4, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_7_p     = [
        # lm.Point(Rbm, theta_p - theta_b / 2 + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rbm, 2*np.pi/p*i + theta_m/2 + theta_g, ms=l3, mode='polar', rotate_instance=True) for i in range(p)
    ]

    magnet_air_slot_8_p    = [
        # lm.Point(Rtm, theta_p - theta_b / 2 + 2 * np.pi / p * i, ms=ks, mode='polar', rotate_instance=True) for i in range(p)
        lm.Point(Rtm, 2*np.pi/p*i + theta_m/2 + theta_g, ms=l4, mode='polar', rotate_instance=True) for i in range(p)
    ]

    # -- Stator Outer Surface
    # stator_outer_surface_p  = [lm.Point(Rout, 0, ms=ks)]
    # [stator_outer_surface_p.extend(lm.rotate([stator_outer_surface_p[0]], angle=[0., 0., (i + 1) * np.pi / 2])) for i in range(2)]

    stator_outer_surface_p = [
        lm.Point(Rout, np.pi * i / 2, ms=l2, mode='polar') for i in range(4)
    ]

    # -- Stator Inner Surface
    shift = theta_t / 2
    # shift = 0
    stator_inner_surface_1_p = [
        # lm.Point(Rsy, theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
        lm.Point(Rsy, theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_2_p = [
        lm.Point(Rsy, theta_sso / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_3_p = [
        lm.Point(Rssi, theta_sso / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_4_p = [
        lm.Point(Rssi, theta_ssi / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_5_p = [
        lm.Point(Rs, theta_ssi / 2 + theta_t * i + shift, ms=l4, mode='polar') for i in range(s)
    ]

    stator_inner_surface_6_p = [
        lm.Point(Rs, theta_t - theta_ssi / 2 + theta_t * i + shift, ms=l4, mode='polar') for i in range(s)
    ]

    stator_inner_surface_7_p = [
        lm.Point(Rssi, theta_t - theta_ssi / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_8_p = [
        lm.Point(Rssi, theta_t - theta_sso / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    stator_inner_surface_9_p = [
        lm.Point(Rsy, theta_t - theta_sso / 2 + theta_t * i + shift, ms=l3, mode='polar') for i in range(s)
    ]

    # -- Domain Boundary
    domain_boundary_p       = [lm.Point(D, 0, ms=l1)]
    for i in range(3):
        domain_boundary_p.extend(lm.rotate([domain_boundary_p[0]], angle = [0., 0., (i + 1) * np.pi / 2]))

    # -- Stator Winding Mid-Section
    stator_winding_mid_p    = [lm.Point(Rssi, 0. + shift, ms=l3)]
    [stator_winding_mid_p.extend(lm.rotate([stator_winding_mid_p[0]], angle = [0., 0., theta_t * (i + 1)])) for i in range(s-1)]

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
    rotor_inner_surface_c  = [
        lm.Curve(rotor_inner_surface_p[i], rotor_inner_surface_p[i+1], origin, curve_type='arc') for i in range((len(rotor_inner_surface_p) - 1))
    ]
    rotor_inner_surface_c.append(
        lm.Curve(rotor_inner_surface_p[-1], rotor_inner_surface_p[0], origin, curve_type='arc')
    )

    # -- Rotor Outer Surface
    rotor_outer_surface_c   = []
    # for i in range(p):
    #     rotor_outer_surface_c.append(
    #         lm.Curve(rotor_outer_surface_p[-1 - i], rotor_outer_surface_p[-2 - i], origin, curve_type='arc')
    #     )

    for i in range(len(rotor_outer_surface_p) - 1):
        rotor_outer_surface_c.append(
            lm.Curve(rotor_outer_surface_p[i], rotor_outer_surface_p[i + 1], origin, curve_type='arc')
        )
    rotor_outer_surface_c.append(lm.Curve(rotor_outer_surface_p[-1], rotor_outer_surface_p[0], origin, curve_type='arc'))



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
            # lm.Curve(magnet_air_slot_3_p[i], magnet_air_slot_4_p[i], physical_group=(2*i+1, 'Magnet Side R ' + str(i+1))),
            lm.Curve(magnet_air_slot_3_p[i], magnet_air_slot_4_p[i]),
            lm.Curve(magnet_air_slot_4_p[i], magnet_air_slot_6_p[i], origin, curve_type='arc'),
            # lm.Curve(magnet_air_slot_6_p[i], magnet_air_slot_5_p[i], physical_group=(2*i+2, 'Magnet Side L ' + str(i+1))),
            lm.Curve(magnet_air_slot_6_p[i], magnet_air_slot_5_p[i]),
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
        lm.Curve(stator_outer_surface_p[i], stator_outer_surface_p[i+1], origin, curve_type='arc') for i in range((len(stator_outer_surface_p) - 1))
    ]
    outer_stator_surface_c.append(
        lm.Curve(stator_outer_surface_p[-1], stator_outer_surface_p[0], origin, curve_type='arc')
    )

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
    ] for i in range(s -  1)])

    stator_core_boundary_c.append([
        lm.Curve(stator_inner_surface_1_p[-1], stator_inner_surface_2_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_2_p[-1], stator_inner_surface_3_p[-1]),
        lm.Curve(stator_inner_surface_3_p[-1], stator_inner_surface_4_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_4_p[-1], stator_inner_surface_5_p[-1]),
        lm.Curve(stator_inner_surface_5_p[-1], stator_inner_surface_6_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_6_p[-1], stator_inner_surface_7_p[-1]),
        lm.Curve(stator_inner_surface_7_p[-1], stator_inner_surface_8_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_8_p[-1], stator_inner_surface_9_p[-1]),
        lm.Curve(stator_inner_surface_9_p[-1], stator_inner_surface_1_p[0], origin, curve_type='arc'),
    ])
    print(len(stator_core_boundary_c[0]))
    stator_core_boundary_c  = list(
        np.array(stator_core_boundary_c).reshape(s*len(stator_core_boundary_c[0]))
    )
    # stator_core_boundary_c.extend([
    #     lm.Curve(stator_inner_surface_1_p[-1], stator_outer_surface_p[-1]),
    #     lm.Curve(stator_outer_surface_p[1], stator_outer_surface_p[2], origin, curve_type='arc'),
    #     lm.Curve(stator_outer_surface_p[0], stator_outer_surface_p[1], origin, curve_type='arc'),
    #     lm.Curve(stator_outer_surface_p[0], stator_inner_surface_1_p[0])
    # ])

    # -- Stator Winding Curves
    stator_windings_c = []
    stator_windings_c.append([
        # lm.Curve(stator_winding_mid_p[0], stator_inner_surface_7_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_7_p[-1], stator_inner_surface_8_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_8_p[-1], stator_inner_surface_9_p[-1]),
        lm.Curve(stator_inner_surface_9_p[-1], stator_inner_surface_1_p[0], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_1_p[0], stator_inner_surface_2_p[0], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_2_p[0], stator_inner_surface_3_p[0]),
        lm.Curve(stator_inner_surface_3_p[0], stator_inner_surface_4_p[0], origin, curve_type='arc'),
        # lm.Curve(stator_inner_surface_4_p[0], stator_winding_mid_p[0], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_4_p[0], stator_inner_surface_7_p[-1], origin, curve_type='arc'),
    ])
    stator_windings_c.extend([[
        # lm.Curve(stator_winding_mid_p[i+1], stator_inner_surface_7_p[i], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_7_p[i], stator_inner_surface_8_p[i], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_8_p[i], stator_inner_surface_9_p[i]),
        lm.Curve(stator_inner_surface_9_p[i], stator_inner_surface_1_p[i+1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_1_p[i+1], stator_inner_surface_2_p[i+1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_2_p[i+1], stator_inner_surface_3_p[i+1]),
        lm.Curve(stator_inner_surface_3_p[i+1], stator_inner_surface_4_p[i+1], origin, curve_type='arc'),
        # lm.Curve(stator_inner_surface_4_p[i+1], stator_winding_mid_p[i+1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_4_p[i+1], stator_inner_surface_7_p[i], origin, curve_type='arc'),
    ] for i in range(s-1)])

    stator_windings_c =  list(
        np.array(stator_windings_c).reshape((s * len(stator_windings_c[0])))
    )

    # -- Air Gap Curves
    air_gap_curves = []
    air_gap_curves.extend([[
        # lm.Curve(stator_winding_mid_p[i], stator_inner_surface_4_p[i], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_4_p[i], stator_inner_surface_5_p[i]),
        lm.Curve(stator_inner_surface_5_p[i], stator_inner_surface_6_p[i], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_6_p[i], stator_inner_surface_7_p[i]),
        # lm.Curve(stator_inner_surface_7_p[i], stator_winding_mid_p[i + 1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_7_p[i], stator_inner_surface_4_p[i + 1], origin, curve_type='arc'),
    ] for i in range(s - 1)])

    air_gap_curves.append([
        # lm.Curve(stator_winding_mid_p[-1], stator_inner_surface_4_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_4_p[-1], stator_inner_surface_5_p[-1]),
        lm.Curve(stator_inner_surface_5_p[-1], stator_inner_surface_6_p[-1], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_6_p[-1], stator_inner_surface_7_p[-1]),
        # lm.Curve(stator_inner_surface_7_p[-1], stator_winding_mid_p[0], origin, curve_type='arc'),
        lm.Curve(stator_inner_surface_7_p[-1], stator_inner_surface_4_p[0], origin, curve_type='arc'),
    ])

    air_gap_curves = list(
        np.array(air_gap_curves).reshape((s * len(air_gap_curves[0])))
    )

    # air_gap_curves.extend(stator_entities_boundary_c[:45])
    # air_gap_curves.append(lm.Curve(stator_winding_mid_p[-1],rotor_outer_surface_p[-1]))
    # air_gap_curves.extend([rotor_outer_surface_c[i] for i in range(8)])
    # air_gap_curves.append(lm.Curve(rotor_outer_surface_p[0], stator_winding_mid_p[0]))

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

    stator_core_outer_surface     = lm.Surface(outer_stator_surface_c)
    # m.add_entity(stator_core_outer_surface)

    stator_core_inner_surface     = lm.Surface(stator_core_boundary_c)
    # m.add_entity(stator_core_inner_surface)

    air_gap_surface = lm.Surface(air_gap_curves)
    # m.add_entity(air_gap_surface)

    rotor_core_surface = []
    rotor_core_surface.extend(rotor_outer_surface_c)

    rotor_core_surface      = lm.Surface(rotor_core_surface)
    # m.add_entity(rotor_core_surface)

    # Magnets and corresponding air slots
    magnet_s, left_air_slot_s, right_air_slot_s = [], [], []
    for i in range(p):
        if i%2 == 0:
            magnet_color = red
        else:
            magnet_color = blue
        magnet_s.append(lm.Surface(magnet_c[4 * i:4 * i + 4], physical_group=(5 + 2*p + i, 'Magnet '+str(i+1)), color=magnet_color))
        left_air_slot_s.append(lm.Surface(left_air_slot_c[4 * i:4 * i + 4], physical_group=(5 + 2*i, 'Left Air Slot '+str(i+1)), color=white))
        right_air_slot_s.append(lm.Surface(right_air_slot_c[4 * i:4 * i + 4], physical_group=(5 + 2*i + 1, 'Right Air Slot '+str(i+1)), color=white))
        # m.add_entity(magnet_s[i])
        # m.add_entity(left_air_slot_s[i])
        # m.add_entity(right_air_slot_s[i])

    rotor_inner_surface_s   = lm.Surface(
        rotor_inner_surface_c,
        physical_group=(3, "Shaft"), 
        color=white
    )
    # m.add_entity(rotor_inner_surface_s)

    # ---------------------------- BOOLEAN SURFACES ----------------------------


    rotor_core_subtract = []
    rotor_core_subtract.extend(magnet_s)
    rotor_core_subtract.extend(left_air_slot_s)
    rotor_core_subtract.extend(right_air_slot_s)
    rotor_core_subtract.append(rotor_inner_surface_s)

    # stator_cut = []
    # stator_cut.extend(right_winding_surfaces)
    # stator_cut.extend(left_winding_surfaces)
    # stator_cut.append(air_gap_surface)
    # stator_cut.extend()

    # stator_core = lm.BooleanSurface(
    #     [stator_core_outer_surface], stator_cut, removeTool=True, operation='subtract',
    # )
    # m.add_entity(stator_core)

    air_gap = lm.BooleanSurface(
        [air_gap_surface], [rotor_core_surface], removeTool=True, operation='subtract', physical_group=(4, 'Air Gap'), color=white
    )
    m.add_entity(air_gap)

    rotor_core = lm.BooleanSurface(
        [rotor_core_surface], rotor_core_subtract, removeTool=False, operation='subtract', physical_group=(1, 'Rotor Core'), color=black
    )
    m.add_entity(rotor_core)

    winding_surfaces = []
    for i in range(s):
        if i%3 == 0: # PHASE B
            winding_color = yellow
            phase = 'B '
        elif i%3 == 1: # PHASE A
            winding_color = green
            phase = 'A '
        elif i%3 == 2: # PHASE C
            winding_color = pink
            phase = 'C '
        winding_surfaces.append(lm.Surface(stator_windings_c[7 * i: 7 * i + 7], physical_group=(5 + 3*p + i, 'Winding ' + phase + str(i+1)), color=winding_color))
        m.add_entity(winding_surfaces[i])

    asdf  = lm.Surface([outer_stator_surface_c],[stator_core_boundary_c], input_type='curve_loops', physical_group=(2, 'Stator Core'), color=black)
    m.add_entity(asdf)

    '''  
    ========= PHYSICAL GROUP ORDERING =========
    Rotor Core (Neodymium):                     1
    Stator Core (Neodymium):                    2
    Shaft (Titanium/Air):                       3
    Air Gap:                                    4
    Air Slots (on the right, n = radially out): 5 + 2*i,            i = 0:(p-1) (these hold odd indices)
    Air Slots (on the left, n = radially out):  5 + 2*i + 1,        i = 0:(p-1) (these hold even indices)
    Magnets:                                    5 + 2*p + i,        i = 0:(p-1)
    Right Stator Windings:                      5 + 3*p + 2*i,      i = 0:(s-1)
    Left Stator Windings:                       5 + 3*p + 2*i  + 1  i = 0:(s-1)
    '''
    # asdf = []
    # asdf.append(air_gap)
    # asdf.append(rotor_core)
    # asdf.extend(stator_cut)

    # stator_core  = lm.BooleanSurface(
    #     [stator_core_outer_surface], asdf, removeTool=True, operation='subtract',
    # )
    # m.add_entity(stator_core)








    # rotor_core_subtract = []
    # rotor_core_subtract.append(rotor_inner_surface_s)
    # rotor_core_subtract.extend(magnet_s)
    # rotor_core_subtract.extend(left_air_slot_s)
    # rotor_core_subtract.extend(right_air_slot_s)
    # rotor_core              = lm.BooleanSurface(
    #     [rotor_core_surface], rotor_core_subtract, removeTool=False, operation='subtract', physical_group=(1, 'Rotor Core')
    # )
    # rotor_core = rotor_core_surface
    # m.add_entity(rotor_core)



    # domain_boundary_surface     = lm.Surface(domain_boundary_c, physical_group = (47, 'Outer Domain'))
    # m.add_entity(domain_boundary_surface)

    # air_gap_surface = lm.Surface(air_gap_curves, physical_group = (46, 'Air Gap'))
    # m.add_entity(air_gap_surface)

    # m.add_all_entities_to_physical_group(geometry_type='curves')



    # FFD PORTION

    mesh_points  = m.get_coordinates(coord_sys='cartesian')
    # print(mesh_points)

    # plt.plot(
    #     [point[0] for point in mesh_points],
    #     [point[1] for point in mesh_points],
    #     'k*'
    # )
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

    # m.add_face(inner_rotor_f[0])

    aaa = lm.Face([magnet_air_slot_1_p[0], magnet_air_slot_2_p[0]], input_type='polar')
    aaa.add_shape_parameter(
        'magnet_point', 
        'r', 
        'constant'
    )
    # m.add_face(aaa)

    bbb = lm.Face([magnet_air_slot_1_p[1], magnet_air_slot_2_p[1]], input_type='polar')
    bbb.add_shape_parameter(
        'magnet_point', 
        'theta', 
        'linear'
    )
    m.add_face(bbb)

    pi_test = lm.Face([magnet_air_slot_3_p[2], magnet_air_slot_6_p[2]], input_type='polar')
    pi_test.add_shape_parameter(
        'inner_magnet_radius',
        'r',
        'constant'
    )
    # m.add_face(pi_test)

    # pi_test = lm.Face([magnet_air_slot_3_p[1]], input_type='polar')
    # pi_test.add_shape_parameter(
    #     'inner_magnet_radius',
    #     'r',
    #     'constant'
    # )
    # m.add_face(pi_test)

    # BELOW TEST FFD FACES FOR NON-ROTATING INSTANCES, LIKE IN THE STATOR
    stator_ffd = lm.Face([stator_inner_surface_1_p[0], stator_inner_surface_2_p[0]], input_type='polar')
    stator_ffd.add_shape_parameter(
        'stator_test_ffd',
        'theta',
        'constant',
    )
    # m.add_face(stator_ffd)

    Ru_test_def_ffd = lm.Face(
        [magnet_air_slot_1_p[0],
            magnet_air_slot_2_p[0],
            magnet_air_slot_3_p[0],
            magnet_air_slot_4_p[0],
            magnet_air_slot_5_p[0],
            magnet_air_slot_6_p[0],
            magnet_air_slot_7_p[0],
            magnet_air_slot_8_p[0]],
        input_type='polar',
    )
    Ru_test_def_ffd.add_shape_parameter(
        'Ru_test_def_ffd',
        'r',
        'constant',
    )
    # m.add_face(Ru_test_def_ffd)

    m.add_all_entities_to_physical_group('curves')

    m.get_node_indices(stator_inner_surface_4_p, print_to_file='A_z_stator_indices') # Getting node indices out for FLUX LINKAGE

    m.assemble(coordinate_system='polar')

    def getInitialEdgeCoords():
        old_edge_coords = m.get_ffd_edge_old_coords(output_type='cartesian')
        # Ru: trim out the origin (x=0,y=0) where there's no nearby (dist<1e-10) nodes in the mesh
        return old_edge_coords[:-2]
    def generateMeshMovement(angle=0.):
        # delta = np.zeros((4 * vars(m)['num_ffd_faces'], 2)) # array of deltas applied to FFD FACES
        # delta[:8, 1] = 0
        # for i in range(4):
        #     delta[2 * i, 0] = angle
        #     delta[2 * i + 1, 0] = -angle
        # #delta[:8, 0] = np.pi/6

        # delta[8:, 1] = 0
        # delta[8:, 0] = 0
        delta = np.zeros((4 * vars(m)['num_ffd_faces'], 2))
        # delta = np.zeros((8, 2))
        for i in range(4):
            delta[i, 1] = -.02 # shifting first magnet + air gaps radially in by 0.02 m
        edge_deltas= m.test_ffd_edge_parametrization_polar(delta,   
                                                    output_type='cartesian')
        # print(edge_deltas)
        # exit()
        return edge_deltas[:-2]

    old_edge_coords = getInitialEdgeCoords()
    edge_deltas     = generateMeshMovement()

    

    init_edge_coords = 'init_edge_coords.txt'
    # f1  = open(init_edge_coords, 'w')
    # for i in range(old_edge_coords.shape[0]):
    #     f.write(old_edge_coords[i] + '\n')
    # f1.close()
    np.savetxt(init_edge_coords, old_edge_coords)

    edge_coord_deltas = 'edge_coord_deltas.txt'
    # f2  = open(edge_coord_deltas, 'w')
    # for i in range(edge_deltas.shape[0]):
    #     f.write(edge_deltas[i] + '\n')
    # f2.close()
    np.savetxt(edge_coord_deltas, edge_deltas)

    return m

# os.system('python3 msh2xdmf.py -d 2 ' + file_name + '.msh')
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


if __name__ == '__main__':
    mesh_object = MotorMeshGenerator(
        rotation_angles=[0], 
        file_name='mesh_files/motor_mesh_new',
        poles=12,
    )   

