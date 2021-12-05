import lsdo_mesh.geometry_classes as lm
import numpy as np

filename = 'magnet_disk_mesh'
m = lm.Mesh(name=filename, popup=True)

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

f3 = lm.Face([magnet_points[1][0], magnet_points[1][1]], input_type='polar')
f3.add_shape_parameter(
    name='magnet_inner_radius_1',
    axis='r',
    def_type = 'constant'
)
m.add_face(f3)

f4 = lm.Face([magnet_points[1][2], magnet_points[1][3]], input_type='polar')
f4.add_shape_parameter(
    name='magnet_inner_radius_1',
    axis='r',
    def_type = 'constant'
)
m.add_face(f4)

m.assemble(coordinate_system='polar')

#print(vars(m).keys())

# import os
# os.system('python3 msh2xdmf.py -d 2 ' + filename + '.msh')

''' ---------------------------------- FFD TEST ---------------------------------- '''

delta = np.zeros((4 * vars(m)['num_ffd_faces'], 2)) # array of deltas applied to FFD FACES
delta[:8, 1] = 0.25
for i in range(4):
    delta[2 * i, 0] = np.pi/36
    delta[2 * i + 1, 0] = -np.pi/36
#delta[:8, 0] = np.pi/6

delta[8:, 1] = 0.25
delta[8:, 0] = 0

# above entries (0:8) correspond to magnet in Q1 shifting by pi/8 ccw and radially out by 0.25
# above entries (8:16) correspond to magnet in Q2 shifting radially out by 0.25


# TODO: maybe we no longer need the `edge_indices` from `lsdo_mesh` as we now use KDTree
# to extract the node indices from FEniCS; and we should only need either `new_edge_coords` 
# or `edge_deltas`, given that we already have `old_edge_coords`.
def getInitialEdgeCoords():
    old_edge_coords = m.get_ffd_edge_old_coords(output_type='cartesian')
    return old_edge_coords
def generateMeshMovement():
    edge_deltas= m.test_ffd_edge_parametrization_polar(delta,   
                                                output_type='cartesian')
    return edge_deltas
# outputs of line 132 are the new edge coordinates, old edge coordinates, the delta array and edge indices

#print(old_edge_coords)
#print(new_edge_coords)
#print(edge_deltas)
#print(generateMeshMovement()[3])

# ------------------ MANUAL TEST ------------------
if False:  # (can ignore this for now)
    ffd_sparse_mat = vars(m)['ffd_face_sps_mat']
    ordered_coords = vars(m)['gmsh_order_point_coords_polar'][:, :2]
    ordered_coords.reshape((len(ordered_coords)*2,))
#    print(ffd_sparse_mat.shape, delta.shape, ordered_coords.shape)
    new_coords = ffd_sparse_mat.dot(delta.reshape((2 * delta.shape[0],))) + ordered_coords.reshape((2 * ordered_coords.shape[0],))
#    print('---')
#    print(ordered_coords)
#    print(new_coords)

#    print(- ordered_coords.reshape((2 * ordered_coords.shape[0],)) + new_coords)

    # EDGE CHECK
    edge_sparse_mat = vars(m)['edge_param_sps_mat']
    # print(edge_sparse_mat.shape)
    new_edge_coords = edge_sparse_mat.dot(new_coords)
    old_edge_coords = edge_sparse_mat.dot(ordered_coords.reshape((2 * ordered_coords.shape[0],)))

    edge_deltas = new_edge_coords - old_edge_coords
