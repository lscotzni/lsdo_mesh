import numpy as np
import csdl
from csdl_om import Simulator

from motor_mesh import MotorMeshGenerator

rotor_rotations     = np.array([
    i * 45 * np.pi / 180 for i in range(2)
])

mesh_objects = []
base_file_name = 'motor_mesh'
mesh_object = MotorMeshGenerator(rotor_rotations, base_file_name)

csdl_mesh_model = vars(mesh_object)['mesh_model']
print(' ------------------------- CSDL MODEL VARIABLE CHECK -------------------------')
# for key in list(vars(csdl_mesh_model).keys()):
#     print('-----')
#     print(key)
#     print(vars(csdl_mesh_model)[key])

sim = Simulator(csdl_mesh_model)
sim['shape_parameter_vec'] = np.zeros((3,))
sim.run()
print(sim['output'])
result = sim['output']
mesh1_length = int(len(result)/2)

asdf  = result[:mesh1_length]  - result[mesh1_length:]
print(asdf)


# print(list(vars(csdl_mesh_model)))
# print(list(vars(csdl_mesh_model).keys()))
exit()

# TESTING MESH PARAMETRIZATION DIFFERENCES
# won't work anymore because of changes
mesh_1 = vars(mesh_objects[0])
mesh_2 = vars(mesh_objects[1])

edge_sparse_mat = mesh_1['edge_param_sps_mat']
num_pts = mesh_1['mesh_nodes'].shape[0]

mesh_1_polar_coords = mesh_1['gmsh_order_point_coords_polar']
mesh_2_polar_coords = mesh_2['gmsh_order_point_coords_polar']

mesh_1_node_coord_test = edge_sparse_mat.dot(
    mesh_1_polar_coords[:,:2].reshape((num_pts * 2))
)

mesh_2_node_coord_test = edge_sparse_mat.dot(
    mesh_2_polar_coords[:,:2].reshape((num_pts * 2))
)

error = mesh_1_node_coord_test - mesh_2_node_coord_test
error_norm = np.linalg.norm(error)
