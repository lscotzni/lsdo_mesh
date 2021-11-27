import lsdo_mesh.geometry_classes as lm
import numpy as np

filename = 'magnet_disk_mesh'
m = lm.Mesh(name=filename, popup=False)

# PARAMETERS
rd = 6.e0 # disk radius
rm_i = 3.5e0 # inner radius of magnet
rm_o = 4.e0 # outer radius of magnet
tm = np.pi / 8 # magnet angular sweep
ti = np.pi / 4
# magnets at 45 degrees + n * pi/2
# disk material: electrical/laminated steel
# magnet material: Neodymium magnet
num_magnets = 4

ks = 0.3

# ----------------- ADDING POINTS -----------------
# ORIGIN
origin = lm.Point(0.0, 0.0, ms=ks)

# DISK BOUNDARY
disk_bound_points = [
    lm.Point(rd, np.pi / 2 * i, ms=ks, mode='polar') for i in range(4)
]
 
# MAGNET POINTS:
magnet_points = [[
    lm.Point(rm_o, ti - tm / 2, ms=ks, mode='polar'),
    lm.Point(rm_o, ti + tm / 2, ms=ks, mode='polar'),
    lm.Point(rm_i, ti + tm / 2, ms=ks, mode='polar'),
    lm.Point(rm_i, ti - tm / 2, ms=ks, mode='polar'),
]]
for i in range(3):
    magnet_points.append(
        lm.rotate(magnet_points[0], angle=[0., 0., np.pi / 2 * (i+1)])
    )

# ----------------- ADDING CURVES -----------------
disk_bound_curves = []
for i in range(len(disk_bound_points)):
    if i == len(disk_bound_points) - 1:
        disk_bound_curves.append(
            lm.Curve(disk_bound_points[i], disk_bound_points[0], origin, curve_type='arc')
        )
    else:
        disk_bound_curves.append(
            lm.Curve(disk_bound_points[i], disk_bound_points[i+1], origin, curve_type='arc')
        )

magnet_curves = []
for i in range(len(magnet_points)):
    magnet_curves.append([
        lm.Curve(magnet_points[i][0], magnet_points[i][1], origin, curve_type='arc'),
        lm.Curve(magnet_points[i][1], magnet_points[i][2]),
        lm.Curve(magnet_points[i][2], magnet_points[i][3], origin, curve_type='arc'),
        lm.Curve(magnet_points[i][3], magnet_points[i][0]),
    ])

# ----------------- ADDING SURFACES -----------------
disk_bound_surface = lm.Surface(disk_bound_curves)

magnet_surfaces = [
    lm.Surface(magnet_curves[i], physical_group=(2+i, 'Magnet ' + str(i+1))) for i in range(len(magnet_curves))
]

disk_surface = lm.BooleanSurface(
    [disk_bound_surface], 
    magnet_surfaces, 
    operation='subtract', 
    removeTool=False,
    physical_group=(1, 'Disk Surface')
)
m.add_entity(disk_surface)

m.add_all_entities_to_physical_group(geometry_type='curves')




''' ---------------------------------- FFD ---------------------------------- '''
f1 = lm.Face([magnet_points[0][0], magnet_points[0][1]], input_type='polar')
f1.add_shape_parameter(
    name='magnet_outer_radius_1',
    axis='r',
    def_type = 'constant'
)
m.add_face(f1)

f2 = lm.Face([magnet_points[0][2], magnet_points[0][3]], input_type='polar')
f2.add_shape_parameter(
    name='magnet_inner_radius_1',
    axis='r',
    def_type = 'linear'
)
m.add_face(f2)

# f3 = lm.Face([magnet_points[0][1]], input_type='polar')
# f3.add_shape_parameter(
#     name='magnet_inner_radius_1',
#     axis='r',
#     def_type = 'constant'
# )
# m.add_face(f3)

# f4 = lm.Face([magnet_points[0][0]], input_type='polar')
# f4.add_shape_parameter(
#     name='magnet_inner_radius_1',
#     axis='r',
#     def_type = 'constant'
# )
# m.add_face(f4)

m.assemble(coordinate_system='polar')
# print(ti - tm / 2)
# print(ti + tm / 2)

# print(vars(f1))

# import os
# os.system('python3 msh2xdmf.py -d 2 ' + filename + '.msh')



