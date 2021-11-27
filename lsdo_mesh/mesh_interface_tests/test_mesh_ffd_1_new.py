import lsdo_mesh.geometry_classes as lm
import numpy as np
import matplotlib.pyplot as plt

m = lm.Mesh(name='test_mesh_ffd_1', popup=True)

ri = 5.
ro = 7.
ts = np.pi / 8.
t0 = np.pi / 4.
ms = 0.5

origin = lm.Point(0., 0., ms=ms)

# Defining circle sector in POLAR (means it should create arcs in cartesian grid)
p1 = lm.Point(ro+1., t0, ms=ms)
p2 = lm.Point(ro+1., t0 + ts, ms=ms)
p3 = lm.Point(ri+1., t0 + ts, ms=ms)
p4 = lm.Point(ri+1., t0, ms=ms)

c1 = lm.Curve(p1, p2, origin, coord_sys='polar')
c2 = lm.Curve(p2, p3, coord_sys='polar')
c3 = lm.Curve(p3, p4, origin, coord_sys='polar')
c4 = lm.Curve(p4, p1, coord_sys='polar')

s1 = lm.Surface([c1, c2, c3, c4], input_type='curves')
# s1_p = lm.Surface([p1, p2, p3, p4], input_type='polygon')
m.add_entity(s1)

# Defining circle sector in CARTESIAN
p5 = lm.Point(ro * np.cos(t0), ro * np.sin(t0), ms=ms)
p6 = lm.Point(ro * np.cos(t0 + ts), ro * np.sin(t0 + ts), ms=ms)
p7 = lm.Point(ri * np.cos(t0 + ts), ri * np.sin(t0 + ts), ms=ms)
p8 = lm.Point(ri * np.cos(t0), ri * np.sin(t0), ms=ms)

c5 = lm.Curve(p5, p6, origin, curve_type='arc')
c6 = lm.Curve(p6, p7)
c7 = lm.Curve(p7, p8, origin, curve_type='arc')
c8 = lm.Curve(p8, p5)

s2 = lm.Surface([c5, c6, c7, c8], input_type='curves')
m.add_entity(s2)

m.assemble()

exit()
# need to figure out how to incorporate the differences in polar vs 
# rectilinear in rest of the objects; the above curves should create
# something like a circle sector in x-y space with the polar argument
# ----------------------------------- FFD -----------------------------------
dr, dt = 0.1, np.pi/20
v1 = lm.Vertex(ro + dr, t0 - dt, mode='polar')
v2 = lm.Vertex(ro + dr, t0 + dt, mode='polar')
v3 = lm.Vertex(ro - dr, t0 + dt, mode='polar')
v4 = lm.Vertex(ro - dr, t0 - dt, mode='polar')

# f1 = lm.Face(p1, input_type='polar')
# f.add_ffd_entity(f1)
# add ffd_entities to Mesh & remove class FFD
# print(vars(f1))
# exit()
# f2 = lm.Face([v1, v2, v3, v4], input_type='polygon')
# f.add_ffd_entity(f2)

# asdf = f.add_shape_parameter('f1_dr', dv1_bounds=[-0.5, 0.5], dv_type='polar')

# make add_dv a method to the ffd Face
# print(asdf)
# exit()
# m.assemble()




# c1 = lm.Curve(..., coord_sys='polar')
# c2 = lm.Curve(..., coord_sys='polar')
# c3 = lm.Curve(..., coord_sys='rectilinear')
# c4 = lm.Curve(..., coord_sys='polar')

# s1 = lm.Surface([c1, c2, c3, c4])
# m.add_entity(s1)