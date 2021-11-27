from lsdo_mesh import geometry_classes as lm
import numpy as np


square_edge     = 10.
ci              = 3.
co              = 7.
instances       = 9
sweep           = np.pi / 20
ms              = 0.25

m               = lm.Mesh(name = 'test_boolean_rotate', popup = True)

# POINTS
origin      = lm.Point(0., 0., ms=ms)

# -- Square
# p1_s        = lm.Point(square_edge, square_edge, ms=ms)
# p2_s        = lm.Point(-square_edge, square_edge, ms=ms)
# p3_s        = lm.Point(-square_edge, 0., ms=ms)
# p4_s        = lm.Point(square_edge, 0., ms=ms)

p1_s        = lm.Point(square_edge, 0, ms=ms)
p2_s        = lm.Point(square_edge, np.pi/2, ms=ms, mode='polar')
p3_s        = lm.Point(square_edge, np.pi, ms=ms, mode='polar')

# -- Slices
p1_a        = lm.Point(ci, 0., ms=ms, mode='polar')
p2_a        = lm.Point(co, 0., ms=ms, mode='polar')
p3_a        = lm.Point(co, sweep, ms=ms, mode='polar')
p4_a        = lm.Point(ci, sweep, ms=ms, mode='polar')

# -- Manual Rotation --
# p1_array    = []
# p2_array    = []
# p3_array    = []
# p4_array    = []

# for i in range(instances):
#     p1_array.append(lm.Point(ci, i * np.pi / instances, ms=ms, mode='polar'))
#     p2_array.append(lm.Point(co, i * np.pi / instances, ms=ms, mode='polar'))
#     p3_array.append(lm.Point(co, i * np.pi / instances + sweep, ms=ms, mode='polar'))
#     p4_array.append(lm.Point(ci, i * np.pi / instances + sweep, ms=ms, mode='polar'))

# -- Using lm.rotate --
p1_array    = [p1_a]
p2_array    = [p2_a]
p3_array    = [p3_a]
p4_array    = [p4_a]

for i in range(instances-1):
    p1_array.extend(lm.rotate([p1_a], angle=[0., 0., (i + 1) * np.pi / instances]))
    p2_array.extend(lm.rotate([p2_a], angle=[0., 0., (i + 1) * np.pi / instances]))
    p3_array.extend(lm.rotate([p3_a], angle=[0., 0., (i + 1) * np.pi / instances]))
    p4_array.extend(lm.rotate([p4_a], angle=[0., 0., (i + 1) * np.pi / instances]))


# CURVES
# -- Square
ct_s        = lm.Curve(p1_s, p2_s, origin, curve_type='arc')
cl_s        = lm.Curve(p2_s,p3_s, origin, curve_type='arc')
cb_s        = lm.Curve(p3_s,origin)
cr_s        = lm.Curve(origin,p1_s)

# -- Slices
c1_array    = []
c2_array    = []
c3_array    = []
c4_array    = []

for i in range(instances):
    # print(i)
    c1_array.append(lm.Curve(p1_array[i],p2_array[i]))
    c2_array.append(lm.Curve(p2_array[i],p3_array[i],origin,curve_type = 'arc'))
    c3_array.append(lm.Curve(p3_array[i],p4_array[i]))
    c4_array.append(lm.Curve(p4_array[i],p1_array[i],origin,curve_type = 'arc'))

# SURFACES
# -- Square
s1_s        = lm.Surface([ct_s, cl_s, cb_s, cr_s])
# m.add_entity(s1_s)

# -- Slices
s_array     = []
for i in range(instances):
    s_array.append(lm.Surface([
        c1_array[i],
        c2_array[i],
        c3_array[i],
        c4_array[i]
    ]))
    # m.add_entity(s_array[i])

final_surface   = lm.BooleanSurface([s1_s], s_array, operation = 'subtract')
m.add_entity(final_surface)
m.assemble()
