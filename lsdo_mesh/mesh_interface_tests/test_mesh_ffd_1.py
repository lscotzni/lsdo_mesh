import lsdo_mesh.geometry_classes as lm
import numpy as np
import matplotlib.pyplot as plt

m = lm.Mesh(name='test_mesh_ffd_1', popup=False)

ri = 1.5
ro = 1.
ts = np.pi / 8.
t0 = np.pi / 4.
ms = 0.2

origin = lm.Point(0., 0., ms=ms)

p1 = lm.Point(ro, t0, ms=ms, mode='polar')
p2 = lm.Point(ro, t0 + ts, ms=ms, mode='polar')
p3 = lm.Point(ri, t0 + ts, ms=ms, mode='polar')
p4 = lm.Point(ri, t0, ms=ms, mode='polar')

c1 = lm.Curve(p1, p2, origin, curve_type='arc')
c2 = lm.Curve(p2, p3)
c3 = lm.Curve(p3, p4, origin, curve_type='arc')
c4 = lm.Curve(p4, p1)

s = lm.Surface([c1, c2, c3, c4], input_type='curves')
# s = lm.Surface([p1, p2, p3, p4], input_type='polygon') # ERROR APPEARS BC EDGE PARAMETRIZATION ONLY EXISTS IN POLAR
m.add_entity(s)


# ----------------------------------- FFD -----------------------------------
dr, dt = 0.05, np.pi/50
v1 = lm.Vertex(ro + dr, t0 - dt, mode='polar')
v2 = lm.Vertex(ro + dr, t0 + dt, mode='polar')
v3 = lm.Vertex(ro - dr, t0 + dt, mode='polar')
v4 = lm.Vertex(ro - dr, t0 - dt, mode='polar')
# print(v1.return_coordinates())

v = [v1, v2, v3, v4]

f1 = lm.Face([p1], input_type='polar', delta=10)
f1.add_shape_parameter(name='p1_radius', axis='r', def_type='constant')
# change degree to type = constant, linear, ...
m.add_face(f1)


# f2 = lm.Face([p2], input_type='polar')
# f2.add_shape_parameter(name='p2_radius', axis='r', def_type='linear')
# m.add_face(f2)

# f34 = lm.Face([p3, p4], input_type='polar', delta=15)
# f34.add_shape_parameter(name='p34_radius', axis='r', def_type='linear')
# m.add_face(f34)

m.assemble(coordinate_system='polar')

mesh_model = m.mesh_model
print(vars(mesh_model))
# vertices = vars(f1)['face_vertices']
# print(vars(vertices[0])) # the face vertices properties hold the cartesian coordinates
# print(f1.return_coordinates())
# print(f1.return_coordinates('polar'))
# print('children (after assemble): ')
# print(vars(f1)['children'])
# print(vars(f2)['children'])

exit()
if False:
    mesh_points = m.get_coordinates(coord_sys='polar')
    print(mesh_points)
    edge_points = m.get_edge_parametrization_test(mesh_points[0], mesh_points[1], coord_sys='polar')
    print(edge_points)



# ----------------------------------- PLOTTING p1 & FFD FACE ----------------------------------- 
vertices = f1.return_coordinates()

plt.figure(1)
for i in range(len(vertices)):
    if i != len(vertices) - 1:
        plt.plot(
            [vertices[i][0], vertices[i+1][0]],
            [vertices[i][1], vertices[i+1][1]]
        )
    else:
        plt.plot(
            [vertices[i][0], vertices[0][0]],
            [vertices[i][1], vertices[0][1]]
        )
plt.plot(
    [0., p1.return_coordinates()[0]],
    [0., p1.return_coordinates()[1]], 
    'k-*'
)

plt.plot([0., vertices[0][0]],[0., vertices[0][1]], 'k-*')
plt.plot([0., vertices[1][0]],[0., vertices[1][1]], 'k-*')

plt.show()

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





# c1 = lm.Curve(..., coord_sys='polar')
# c2 = lm.Curve(..., coord_sys='polar')
# c3 = lm.Curve(..., coord_sys='rectilinear')
# c4 = lm.Curve(..., coord_sys='polar')

# s1 = lm.Surface([c1, c2, c3, c4])
# m.add_entity(s1)