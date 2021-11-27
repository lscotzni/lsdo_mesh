from lsdo_mesh import geometry_classes as lm
import numpy as np

'''
TWO RECTANGES CONNECTED
'''

x = 5.
y = 7.

m = lm.Mesh('test_connected_rectangles',popup = True)

# POINTS ----------------------------- 
left_bl_p   = lm.Point(0.,0.)
left_tl_p   = lm.Point(0,y)
left_tr_p   = lm.Point(x,y)
left_br_p   = lm.Point(x,0)
print(vars(left_bl_p))

new_point = lm.rotate([left_br_p],angle = [0.,0.,np.pi/2],input_type = 'points')
print(new_point)
# exit()

mid_bl_p   = lm.Point(x,0)
mid_tl_p   = lm.Point(x,y)
mid_tr_p   = lm.Point(2*x,y)
mid_br_p   = lm.Point(2*x,0)

right_bl_p   = lm.Point(2*x,0)
right_tl_p   = lm.Point(2*x,y)
right_tr_p   = lm.Point(3*x,y)
right_br_p   = lm.Point(3*x,0)

rightright_bl_p   = lm.Point(3*x,0)
rightright_tl_p   = lm.Point(3*x,y)
rightright_tr_p   = lm.Point(4*x,y)
rightright_br_p   = lm.Point(4*x,0)

top_bl_p   = lm.Point(x,y)
top_tl_p   = lm.Point(x,2*y)
top_tr_p   = lm.Point(2*x,2*y)
top_br_p   = lm.Point(2*x,y)

bot_bl_p   = lm.Point(x,0)
bot_tl_p   = lm.Point(x,-y)
bot_tr_p   = lm.Point(2*x,-y)
bot_br_p   = lm.Point(2*x,0)

# CURVES ----------------------------- 

left_l_c    = lm.Curve(left_bl_p,left_tl_p)
left_t_c    = lm.Curve(left_tl_p,left_tr_p) 
left_r_c    = lm.Curve(left_tr_p,left_br_p)
left_b_c    = lm.Curve(left_br_p,left_bl_p)

mid_l_c    = lm.Curve(mid_bl_p,mid_tl_p)
mid_t_c    = lm.Curve(mid_tl_p,mid_tr_p)
mid_r_c    = lm.Curve(mid_tr_p,mid_br_p)
mid_b_c    = lm.Curve(mid_br_p,mid_bl_p)

right_l_c    = lm.Curve(right_bl_p,right_tl_p)
right_t_c    = lm.Curve(right_tl_p,right_tr_p)
right_r_c    = lm.Curve(right_tr_p,right_br_p)
right_b_c    = lm.Curve(right_br_p,right_bl_p)

rightright_l_c    = lm.Curve(rightright_bl_p,rightright_tl_p)
rightright_t_c    = lm.Curve(rightright_tl_p,rightright_tr_p)
rightright_r_c    = lm.Curve(rightright_tr_p,rightright_br_p)
rightright_b_c    = lm.Curve(rightright_br_p,rightright_bl_p)

top_l_c    = lm.Curve(top_bl_p,top_tl_p)
top_t_c    = lm.Curve(top_tl_p,top_tr_p)
top_r_c    = lm.Curve(top_tr_p,top_br_p)
top_b_c    = lm.Curve(top_br_p,top_bl_p)

bot_l_c    = lm.Curve(bot_bl_p,bot_tl_p)
bot_t_c    = lm.Curve(bot_tl_p,bot_tr_p)
bot_r_c    = lm.Curve(bot_tr_p,bot_br_p)
bot_b_c    = lm.Curve(bot_br_p,bot_bl_p)

# SURFACES ----------------------------- 

left_square = lm.Surface(
    [left_l_c,left_t_c,left_r_c,left_b_c]
)
m.add_entity(left_square)

mid_square = lm.Surface(
    [mid_l_c, mid_t_c, mid_r_c, mid_b_c]
)
m.add_entity(mid_square)

right_square = lm.Surface(
    [right_l_c, right_t_c, right_r_c, right_b_c]
)
m.add_entity(right_square)

rightright_square = lm.Surface(
    [rightright_l_c, rightright_t_c, rightright_r_c, rightright_b_c]
)
m.add_entity(rightright_square)

top_square = lm.Surface(
    [top_l_c, top_t_c, top_r_c, top_b_c]
)
m.add_entity(top_square)

bot_square = lm.Surface(
    [bot_l_c, bot_t_c, bot_r_c, bot_b_c]
)
m.add_entity(bot_square)

m.assemble()

# print('--------- UNCHANGED RESULTS ---------')
# print(vars(m)['point_coordinates'])
# print(vars(m)['curves'])
# print(vars(m)['curve_indices'])
# print(vars(m)['curve_type'])
# print(vars(m)['surfaces'])
# print(vars(m)['surface_indices'])
# print(vars(m)['surface_type'])