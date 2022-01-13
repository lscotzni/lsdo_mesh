import numpy as np
import csdl
from csdl import Model
from csdl_om import Simulator

# from motor_mesh_new import MotorMeshGenerator
from motor_mesh import MotorMeshGenerator
# from post_processing_dir.efficiency_model import EfficiencyModel
# from post_processing_dir.flux_linkage_model import FluxLinkageModel
# from post_processing_dir.electrical_model import ElectricalModel
# from post_processing_dir.mass_model import MassModel

rotor_rotations     = np.array([
    i * 45 * np.pi / 180 for i in range(2)
])

test = False
if test:
    p = 4
    base_file_name = 'motor_mesh_test'
else:
    p = 12
    base_file_name = 'motor_mesh'

mesh_objects = MotorMeshGenerator(
    rotation_angles=rotor_rotations, 
    file_name=base_file_name,
    poles=p,
) # outputs list of mesh model instances

print(mesh_objects)


exit()
print(' ------------------------- CSDL MODEL VARIABLE CHECK -------------------------')
# for key in list(vars(csdl_mesh_model).keys()):
#     print('-----')
#     print(key)
#     print(vars(csdl_mesh_model)[key])
print(csdl_mesh_model)

sim = Simulator(csdl_mesh_model)
# sim['shape_parameter_vec'] = np.ones((3,))
sim.run()

print(sim['new_mesh_points_1'] - sim['new_mesh_points_2'])
print(sim['new_edge_nodes_1'] - sim['new_edge_nodes_2'])


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
