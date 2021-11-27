from lsdo_mesh import mesh_classes as lm

'''
SIMPLE EXAMPLE OF A CIRCLE IN A SQUARE
'''

m = lm.Mesh("test_square_circle_holes",popup = True) # CALL TO INITALIZE GMSH

square_side         = 6.
circle_radius       = 1.
circle_center       = 1.5

# SQUARE POINTS/LINES
pt_bot_l    = lm.Point(-square_side/2.,-square_side/2.,ms = 0.05)
pt_bot_r    = lm.Point(square_side/2.,-square_side/2.,ms = 0.05)
pt_top_r    = lm.Point(square_side/2.,square_side/2.,ms = 0.05)
pt_top_l    = lm.Point(-square_side/2.,square_side/2.,ms = 0.05)

cu_left     = lm.Curve(pt_bot_l, pt_top_l)
cu_top      = lm.Curve(pt_top_l, pt_top_r)
cu_right    = lm.Curve(pt_top_r, pt_bot_r)
cu_bot      = lm.Curve(pt_bot_r, pt_bot_l)

# cu_left_dup     = lm.Curve(pt_bot_l, pt_top_l)
# cu_top_dup      = lm.Curve(pt_top_l, pt_top_r)
# cu_right_dup    = lm.Curve(pt_top_r, pt_bot_r)
# cu_bot_dup      = lm.Curve(pt_bot_r, pt_bot_l)

# TOP CIRCLE POINTS/LINES
pt_center_t   = lm.Point(0.0,circle_center,ms = 0.05)

pt_top_c_t    = lm.Point(0.,circle_radius + circle_center,ms = 0.05)
pt_left_c_t   = lm.Point(-circle_radius,circle_center,ms = 0.05)
pt_bot_c_t    = lm.Point(0.,-circle_radius + circle_center,ms = 0.05)
pt_right_c_t  = lm.Point(circle_radius,circle_center,ms = 0.05)

cu_top_l_t    = lm.Curve(pt_top_c_t,pt_left_c_t,pt_center_t,curve_type = 'arc')
cu_bot_l_t    = lm.Curve(pt_left_c_t,pt_bot_c_t,pt_center_t,curve_type = 'arc')
cu_bot_r_t    = lm.Curve(pt_bot_c_t,pt_right_c_t,pt_center_t,curve_type = 'arc')
cu_top_r_t    = lm.Curve(pt_right_c_t,pt_top_c_t,pt_center_t,curve_type = 'arc')

# BOTTOM CIRCLE POINTS/LINES
pt_center_b   = lm.Point(0.0,-circle_center,ms = 0.05)

pt_top_c_b    = lm.Point(0.,circle_radius-circle_center,ms = 0.05)
pt_left_c_b   = lm.Point(-circle_radius,-circle_center,ms = 0.05)
pt_bot_c_b    = lm.Point(0.,-circle_radius-circle_center,ms = 0.05)
pt_right_c_b  = lm.Point(circle_radius,-circle_center,ms = 0.05)

cu_top_l_b    = lm.Curve(pt_top_c_b,pt_left_c_b,pt_center_b,curve_type = 'arc')
cu_bot_l_b    = lm.Curve(pt_left_c_b,pt_bot_c_b,pt_center_b,curve_type = 'arc')
cu_bot_r_b    = lm.Curve(pt_bot_c_b,pt_right_c_b,pt_center_b,curve_type = 'arc')
cu_top_r_b    = lm.Curve(pt_right_c_b,pt_top_c_b,pt_center_b,curve_type = 'arc')

# cu_top_l_dup    = lm.Curve(pt_top_c,pt_left_c,pt_center,curve_type = 'arc')
# cu_bot_l_dup    = lm.Curve(pt_left_c,pt_bot_c,pt_center,curve_type = 'arc')
# cu_bot_r_dup    = lm.Curve(pt_bot_c,pt_right_c,pt_center,curve_type = 'arc')
# cu_top_r_dup    = lm.Curve(pt_right_c,pt_top_c,pt_center,curve_type = 'arc')

# SURFACES
su_square_c   = lm.Surface(
    cu_left,cu_top,cu_right,cu_bot,
    input_type = 'curves'
)
# m.add_entity(su_square_c)


# su_square_c_dup   = lm.Surface(
#     cu_left_dup,cu_top_dup,cu_right_dup,cu_bot_dup,
#     input_type = 'curves'
# )
# m.add_entity(su_square_c_dup)


su_circle_t   = lm.Surface(
    cu_top_l_t,cu_bot_l_t,cu_bot_r_t,cu_top_r_t,
    input_type = 'curves'
)

su_circle_b   = lm.Surface(
    cu_top_l_b,cu_bot_l_b,cu_bot_r_b,cu_top_r_b,
    input_type = 'curves'
)


# su_circle_dup   = lm.Surface(
#     cu_top_l_dup,cu_bot_l_dup,cu_bot_r_dup,cu_top_r_dup,
#     input_type = 'curves'
# )
# m.add_entity(su_circle_dup)

# su2 = lm.subtract(su_square_c, su_square)
# print(su_square_c)
# s3 = lm.Boolean([su_square_c],[su_circle_t],operation = 'subtract')
# m.add_entity(s3)

s4 = lm.Boolean([su_square_c],[su_circle_t,su_circle_b],operation = 'subtract')
m.add_entity(s4)

m.assemble()




# # 
# np.dot(a, b)
# a.dot(b)

# # -----------------
# s1 = lm.Surface(c1, c2, c3)
# s2 = lm.Surface(c1, c2, c3, c4)
# s8 = lm.Surface(pt1, pt2, pt3, pt4, type='from_points')
# s3 = lm.Surface(s1, s2, type='subtraction') # one way to do it, but I don't recommend it
# s11, s12, s13 = lm.subtract(s1, s2) # better way - let's go with this one
# s7 = s1.subtract(s2)
# s5 = s1 - s2 # not as clear - how to handle intersections and fragments
# s6 = s1 + s2