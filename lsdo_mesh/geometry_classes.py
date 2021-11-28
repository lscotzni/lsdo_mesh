import gmsh
import numpy as np
import scipy as sp
from scipy.linalg import norm
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix
import sys

from lsdo_mesh.mesh_ffd_objects import *
# from lsdo_mesh.mesh_ffd_objects_new import *
from lsdo_mesh.remove_duplicates import remove_duplicate_points, remove_duplicate_curves

from lsdo_mesh.geometry_operations import rotate

# ------------------------------------------- MESH -------------------------------------------
class Mesh(object):

    def __init__(self, name='Mesh', popup=False):

        self.name       = name
        self.popup      = popup
        
        self.top_entities           = []

        self.point_coordinates      = []
        self.point_mesh_size        = []

        self.curves                 = []
        self.curve_indices          = []
        self.curve_type             = []
        self.curve_coord_sys        = []

        self.surfaces               = []
        self.surface_indices        = []
        self.surface_type           = []

        self.boolean_operations     = []
        self.boolean_entities       = []
        self.boolean_object_indices = []
        self.boolean_tool_indices   = []
        self.boolean_remove_object  = []
        self.boolean_tool_object    = []
        self.boolean_parameters     = []

        self.point_physical_groups              = []
        self.curve_physical_groups              = []
        self.surface_physical_groups            = []    
        self.boolean_surface_physical_groups    = []

        # Flags for if we want to add all entities of a dimension to a physical group
        self.all_points_physical_group = None
        self.all_curves_physical_group = None
        self.all_surfaces_physical_group = None

        self.ffd_faces              = []
        self.ffd_face_control_pts   = []
        self.ffd_sparse_param_mat   = []


    # --------------------- GEOMETRY/MESH ---------------------
    def add_entity(self, entity=None):
        # coordinates is the indicator for either polar or cartesian for full entity
        self.top_entities.append(entity)

    def add_point(self, x, y, z, mesh_size):
        self.point_coordinates.append([x, y, z])
        self.point_mesh_size.append(mesh_size)
        return len(self.point_coordinates) - 1

    def add_curve(self, children, curve_type, coord_sys=None, physical_group=False):
        self.curves.extend([point.id for point in children])
        if self.curve_indices == []: 
            self.curve_indices.append([0,  len(children)])
        else:
            self.curve_indices.append([
                self.curve_indices[-1][-1],
                self.curve_indices[-1][-1] + len(children),
            ])
        # The other way to set up curve_indices without this "if" clause is by
        # taking the length before and after extending curve_indices (may be cleaner)
        
        self.curve_type.append(curve_type)
        self.curve_coord_sys.append(coord_sys)

        if type(physical_group) is not tuple:
            if physical_group is not False:
                raise KeyError('Physical groups must be tuples denoted as (int, string).')
        self.curve_physical_groups.append(physical_group)
        
        return len(self.curve_type) - 1

    def add_surface(self, children, physical_group=False):
        self.surfaces.extend([curve.id for curve in children])
        if self.surface_indices == []:
            self.surface_indices.append([0, len(children)])
        else:
            self.surface_indices.append([
                self.surface_indices[-1][-1],
                self.surface_indices[-1][-1] + len(children),
            ])
        # physical_group holds (int, string) for index and name

        if type(physical_group) is not tuple:
            if physical_group is not False:
                raise KeyError('Physical groups must be tuples denoted as (int, string).')
        self.surface_physical_groups.append(physical_group)

        return len(self.surface_indices) - 1

    def add_boolean(self, children, properties, physical_group=False):
        self.boolean_operations.append(properties[2])
        self.boolean_entities.extend([entity.id for entity in children])
        if self.boolean_object_indices == []:
            self.boolean_object_indices.append([0, properties[0]])
        else:
            self.boolean_object_indices.append([
                self.boolean_tool_indices[-1][-1],
                self.boolean_tool_indices[-1][-1] + properties[0]
            ])
        if self.boolean_tool_indices == []:
            self.boolean_tool_indices.append([properties[0], properties[0] + properties[1]])
        else:
            self.boolean_tool_indices.append([
                self.boolean_object_indices[-1][-1],
                self.boolean_object_indices[-1][-1] + properties[1]
            ])
        self.boolean_remove_object.append(properties[3])
        self.boolean_tool_object.append(properties[4])
        self.boolean_parameters.append(properties[2:])
        
        if type(physical_group) is not tuple:
            if physical_group is not False:
                raise KeyError('Physical groups must be tuples denoted as (int, string).')
        self.boolean_surface_physical_groups.append(physical_group)

        return len(self.boolean_operations) - 1

    def add_physical_group(self, entities):

        if not all([isinstance(entities[0], entity) for entity in entities[1:]]):
            raise KeyError('Inputs must be the same type.')


    def add_all_entities_to_physical_group(self, geometry_type=None):

        if geometry_type is 'points':
            raise KeyError('Not implemented with points yet.')
        elif geometry_type is 'curves':
            dim = 1
            self.all_curves_physical_group = True
        elif geometry_type is 'surfaces':
            raise KeyError('Not implemented with surfaces yet.')

    # --------------------- FFD ---------------------
    def add_face(self, face):
        self.ffd_faces.append(face)

    def assemble_ffd_parametrization(self, coordinate_system='cartesian'):
        
        print(' ============ ASSEMBLING FFD FACE PARAMETRIZATION ============ ')

        # num_points, dim = self.gmsh_order_point_coords.shape[0], self.gmsh_order_point_coords.shape[1]
        num_points, dim = self.gmsh_order_point_coords.shape[0], 2
        self.num_ffd_faces   = len(self.ffd_faces)
        num_ffd_face_coords = 4 * dim # number of coordinate components stored (2D is x1, x2, etc.)'
        print(self.mesh_nodes)
        print(self.gmsh_order_point_coords)
        # exit()

        # FFD PARAMETRIZATION MATRIX ()
        sparse_row, sparse_col, sparse_val = [], [], []
        # self.ffd_face_sps_mat = lil_matrix(
        #     (num_points * dim, num_ffd_face_coords * self.num_ffd_faces)
        # )
        # num col = # of FFD FACES * 4 * 2 
        # P_new = SPS_MAT * (V - V0) + P0
        ffd_face_control_pts = np.zeros((num_ffd_face_coords * self.num_ffd_faces, ))

        for face_ind, face in enumerate(self.ffd_faces):
            # Coordinates of vertices of individual ffd face
            face_coords = np.array(face.return_coordinates(coordinate_system)) # reshape((4 * dim)) to vectorize

            # physical coordinates of ffd face vertices
            # may be good to hold this information in the vertex object
            # i.e. x1 = [x1_min, x1_max]; x2 = [x2_min, x2_max]
            P00 = [np.min(face_coords[:,0]), np.min(face_coords[:,1])]
            P01 = [np.min(face_coords[:,0]), np.max(face_coords[:,1])]
            P10 = [np.max(face_coords[:,0]), np.min(face_coords[:,1])]
            P11 = [np.max(face_coords[:,0]), np.max(face_coords[:,1])]
            
            face_coords = np.array([P00, P10, P01, P11])

            start = num_ffd_face_coords * face_ind
            end = num_ffd_face_coords * face_ind + num_ffd_face_coords
            # in the order of P00, P10, P01, P11
            ffd_face_control_pts[start:end] = face_coords.reshape((dim*len(face_coords),))

            # loop for the number of children within the face
            embedded_points = vars(face)['embedded_points'] 

            for point_ind, point in enumerate(embedded_points):
                point_coords_to_reorder = np.array(point.return_coordinates(output_type='cartesian'))
                print('coord:', point_coords_to_reorder)
                point_coords = point.return_coordinates(coordinate_system)[:2]
                # print('polar coord:', point_coords)
    
                # might be good to keep this comparison in cartesian coordinates
                # polar coordinates are tough b/c np.arctan2() returns angle in range [-pi, pi]
                # pi and -pi represent the same thing in our case but the sign will break the comparison
                index = np.where(
                    np.linalg.norm(
                        point_coords_to_reorder - self.gmsh_order_point_coords,
                        axis=1
                    ) < 1.e-6
                )
                
                [u, v] = [
                    (point_coords[0] - P00[0]) / (P11[0] - P00[0]),
                    (point_coords[1] - P00[1]) / (P11[1] - P00[1]),
                ]

                for i in range(4):
                    sparse_row.extend(
                        np.arange(dim*index[0][0], dim*index[0][0] + dim)
                    )
                sparse_col.extend(np.arange(start, end, dtype=int))
                sparse_val.extend([
                    (1 - u) * (1 - v),  # P00
                    (1 - u) * (1 - v),
                    u * (1 - v),        # P10
                    u * (1 - v),
                    (1 - u) * v,        # P01
                    (1 - u) * v,
                    u * v,              # P11
                    u * v,
                ])

        # self.ffd_face_sps_mat[sparse_row, sparse_col] = sparse_val
        # self.ffd_face_sps_mat.tocsc()
        self.ffd_face_sps_mat = csc_matrix(
            (sparse_val, (sparse_row, sparse_col)),
            shape=(num_points * dim, num_ffd_face_coords * self.num_ffd_faces)
        )

        asdf = self.ffd_face_sps_mat.dot(ffd_face_control_pts) # check for whether FFD faces return original points
        # other ways to do the above dot product:
        # asdf = np.dot(ffd_face_sps_mat, ffd_face_control_pts)
        # asdf = ffd_face_sps_mat @ ffd_face_control_pts 
        print('FFD Parametrization check:')
        print(asdf)
        print(asdf.shape)
        # exit()

            
    def assemble_edge_parametrization(self, coordinate_system='cartesian'):
        print(' ============ ASSEMBLING EDGE PARAMETRIZATION ============ ')

        # SYSTEM LOOKS AS SUCH
        '''
        MATRIX OF SIZE (2 * # OF NODES/DEFINED POINTS X 2 * # OF EDGE NODES)
        '''
        # ======================== QUESTION FOR DR. HWANG / FIXES FOR THE FUTURE: ========================
        # problem with np.arctan2 and polar coordinates is that pi and -pi are virtually the same coordinate
        # distinction has not been made yet here so we need to add a check to the parametrization
        # for curves where the other point has a negative theta value, this means we must use -pi
        # for curves where the other point has a positive theta value, we use pi
        # gmsh will NOT allow curves to be made further than pi, so this ensures consistency of signs
        num_edge_nodes = max(self.edge_node_indices[-1]) 

        # the origin is ALWAYS the last node set according to gmsh so it is not considered here
        # the nodes are ordered from 1 to number of points
        [num_points, dim] = self.mesh_nodes.shape
        num_edges = len((self.edge_node_indices))

        sparse_row, sparse_col, sparse_val = [], [], []
        # self.edge_param_sps_mat = lil_matrix(
        #     (int((num_edge_nodes + 1) * dim), int((num_points) * dim))
        # )
        
        # IDENTITY MAP FOR NODES CORRESPONDING TO USER-DEFINED MESH POINTS
        # KEY NOTE HERE: THE ORIGIN IS ALWAYS THE LAST NODE, SO WE WILL NOT CONSIDER IT
        # THUS, THE SET OF NODES DEFINING EXACTLY THE POINTS IS EQUAL TO 1: num_points - 1, 
        for i in range(1, 2 * (num_points - 1) + 1):
            sparse_row.append(i-1)
            sparse_col.append(i-1)
            sparse_val.append(1)
        
        sparse_row.extend([num_edge_nodes * dim, num_edge_nodes * dim + 1])
        sparse_col.extend([(num_points - 1) * dim, (num_points - 1) * dim + 1])
        sparse_val.extend([1, 1])
        
        # use self.edge_boundary_nodes and self.edge_node_indices
        for i, edge in enumerate(self.edge_node_indices):
            edge_nodes = self.edge_node_coords[i]
            first_edge = edge_nodes[0]

            # DETERMINING IF PARMAETRIZATION IS DONE IN THETA OR R (correspondingly x and y for cartesian)
            # this can ONLY handle rectilinear edges in either polar or cartesian
            if all([np.linalg.norm(node[0] - first_edge[0]) < 1.e-6 for node in edge_nodes[1:]]):
                # CONSTANT THETA CURVE
                var_dim = 1 # THIS MEANS USE RADIUS TO PARAMETRIZE
                
            elif all([np.linalg.norm(edge_nodes[j][0] - edge_nodes[-j-3][0] - np.pi) < 1.e-8 for j in range(int(np.floor((len(edge_nodes) - 2)/2) - 1))]) and abs(first_edge[0] - np.pi) < 1.e-6:
                # CASE OF LINE PASSING THROUGH ORIGIN; CURRENTLY PRODUCES NAN IN PARAMETRIZATION B/C ENDPOINTS HAVE SAME RADIUS
                var_dim = 1

            elif all([np.linalg.norm(node[1] - first_edge[1]) < 1.e-6 for node in edge_nodes[1:]]):
                # CONSTANT RADIUS CURVE
                var_dim = 0  # THIS MEANS USE THETA TO PARAMETRIZE
            
            pi_sign_mismatch = 0
            pi_sign_start, pi_sign_end = 0, 0
            if np.pi in edge_nodes[-2] or np.pi in edge_nodes[-1]:
                [pi_ind, pi_ind_dim] = np.where(abs(np.pi - edge_nodes) < 1e-8)

                # if isinstance(pi_ind, int): # this condition is never true since pi_ind is always a list, but still works (with positive u outside of [0,1])
                if len(pi_ind) == 1:
                    if np.sign(edge_nodes[pi_ind, pi_ind_dim]) != np.sign(edge_nodes[pi_ind - 1, pi_ind_dim]):
                        pi_sign_mismatch = 1
                        if pi_ind[0] == len(edge_nodes) - 1:
                            pi_sign_end = 1
                        elif pi_ind[0] == len(edge_nodes) - 2:
                            pi_sign_start = 1

            start = int(edge[-2] - 1)
            end = int(edge[-1] - 1)
            num_internal_nodes = len(edge) - 2
            # print(start, end, num_internal_nodes)
            u = []
            for j in range(num_internal_nodes):
                u.append(
                    (edge_nodes[j, var_dim] - edge_nodes[-2, var_dim] * (-1)**pi_sign_start) / \
                    (edge_nodes[-1, var_dim] * (-1)**pi_sign_end - edge_nodes[-2, var_dim] * (-1)**pi_sign_start)
                )

                if u[j] < 0 or u[j] > 1:
                    print('parametrization (u) outside of bounds between 0 and 1; u = ', u[j])

                if 1 - u[j] < 0 or 1 - u[j] > 1:
                    print('parametrization converse (1-u) outside of bounds between 0 and 1; 1 - u = ', 1 - u[j])

                sparse_row.extend(
                    [int(2 * (edge[j] - 1))] * 2 + [int(2 * (edge[j] - 1) + 1)] * 2,
                )
                sparse_col.extend([
                    2 * start,
                    2 * end,
                    2 * (start) + 1, 
                    2 * (end) + 1,
                ])
                if pi_sign_mismatch:
                    # this is where we get a negative parametric value (essentially saying u * pi = -u * -pi)
                    sparse_val.extend([
                        (1 - u[j]) * (-1)**pi_sign_start,
                        u[j], 
                        (1 - u[j]) * (-1)**pi_sign_end,
                        u[j],
                    ])
                else:
                    sparse_val.extend([1 - u[j], u[j]] * 2)

        # self.edge_param_sps_mat[sparse_row, sparse_col] = sparse_val
        # self.edge_param_sps_mat.tocsc()
        self.edge_param_sps_mat = csc_matrix(
            (sparse_val, (sparse_row, sparse_col)),
            shape=(int((num_edge_nodes + 1) * dim), int((num_points) * dim)),
        )

        print(self.edge_param_sps_mat.shape)
        print(self.gmsh_order_point_coords_polar[:,:2].shape)
        # TESTING OUTPUT OF APPLYING PARAMETRIZATION MATRIX TO ORIGINAL POINTS
        node_coords_test = self.edge_param_sps_mat.dot(self.gmsh_order_point_coords_polar[:,:2].reshape((int((num_points) * dim),)))
        
        # Checking if the application of projection onto our point coordinates properly returns the nodes
        # along edges by looking at error norm with edge node coordinates from GMSH
        error_norm_array = []
        for i, edge_nodes in enumerate(self.edge_node_indices):
            edge_node_coords = np.array(self.edge_node_coords[i])[:,:2].reshape((2 * len(edge_nodes),))
            for j, node in enumerate(edge_nodes):
                error = np.array(node_coords_test[int(2*(node-1)):int(2*(node-1) + 2)], dtype=float) - \
                    np.array(edge_node_coords[int(2*(j)):int(2*(j) + 2)], dtype=float)
                error_norm_array.append(np.linalg.norm(error)) 
                # error_norm_array.append(norm(error)) # USING NORM FROM SCIPY (FAILS WITH THE nan)

        high_error_ind = np.where(np.abs(error_norm_array) > 1e-8)

        # UNCOMMENT FIRST 3 LINES BELOW TO LOOK AT ERROR
        print('error norm array: ', error_norm_array)
        print('high error norm locations: ', high_error_ind)
        print('norm of error norm array: ', np.linalg.norm(error_norm_array))
        # print(node_coords_test.reshape((int(len(node_coords_test)/2), dim)))

    def assemble_mesh_parametrization(self, coordinate_system='cartesian'):
        pass

    def create_csdl_model(self):
        pass

    # --------------------- ASSEMBLE ---------------------
    def assemble(self, coordinate_system='cartesian'):
        self.assemble_mesh(coordinate_system=coordinate_system)

        self.get_coordinates(coord_sys=coordinate_system) # returns self.mesh_nodes for entire mesh

        self.assemble_ffd_parametrization(coordinate_system=coordinate_system)

        self.assemble_edge_parametrization(coordinate_system=coordinate_system)

        self.assemble_mesh_parametrization(coordinate_system=coordinate_system)

        self.create_csdl_model() # this contains the ffd face/edge & mesh movement csdl model classes
        # spits out the csdl variable containing mesh coordinates
        
    def assemble_mesh(self, coordinate_system='cartesian'):
        # (1) recursively assemble entities data structures
        # (2) remove duplicates
        # (3) gmsh
        # initialize gmsh
        # # create points, curves, surfaces in order (ignoring surfaces defined as a result of boolean operations)
        # # run the boolean operations in order
        # # finalize gmsh 

        print('Starting mesh assembly')

        for entity in self.top_entities:
            entity.assemble(self)

        # converting to numpy arrays
        self.point_coordinates          = np.array(self.point_coordinates)
        self.point_mesh_size            = np.array(self.point_mesh_size)
        self.point_physical_groups      = np.array(self.point_physical_groups)

        self.curves                     = np.array(self.curves)
        self.curve_indices              = np.array(self.curve_indices)
        self.curve_type                 = np.array(self.curve_type)
        self.curve_coord_sys            = np.array(self.curve_coord_sys)
        # self.curve_physical_groups      = np.array(self.curve_physical_groups)

        self.surfaces                   = np.array(self.surfaces)
        self.surface_indices            = np.array(self.surface_indices)
        self.surface_type               = np.array(self.surface_type)
        # self.surface_physical_groups    = np.array(self.surface_physical_groups)

        self.boolean_operations         = np.array(self.boolean_operations)
        self.boolean_entities           = np.array(self.boolean_entities)
        self.boolean_object_indices     = np.array(self.boolean_object_indices)
        self.boolean_tool_indices       = np.array(self.boolean_tool_indices)
        self.boolean_remove_object      = np.array(self.boolean_remove_object)
        self.boolean_tool_object        = np.array(self.boolean_tool_object)
        self.boolean_parameters         = np.array(self.boolean_parameters)
        
        # code block to execute conversion of points into rectilinear space
        # from polar; this also means converting the proper array values
        # BASIC IDEA: turn curve_coord_sys to arrays of zeros
        for i, cs in enumerate(self.curve_coord_sys):
            if cs == 0 or cs is None:
                continue
            elif cs == 1:
                curve_ind_start = self.curve_indices[i,0]
                curve_ind_end = self.curve_indices[i,1]
                for ind in np.arange(curve_ind_start, curve_ind_end):
                    cu = self.curves[ind]
                    r, theta = self.point_coordinates[cu,0], self.point_coordinates[cu,1]
                    self.point_coordinates[cu] = [
                        r * np.cos(theta),
                        r * np.sin(theta),
                        0.0
                    ]
                self.curve_coord_sys[i] = 0

        # removing duplicate points
        print('Starting removal of duplicate points.')
        self.point_coordinates, self.point_mesh_size, new_curves_temp = remove_duplicate_points(self.point_coordinates,self.point_mesh_size,self.curves)
        print('Completed removal of duplicate points.')

        # removing duplicate curves
        print('Starting removal of duplicate curves.')
        self.curves, self.curve_indices, self.curve_type, self.curve_physical_groups, self.surfaces, unique_surfaces = remove_duplicate_curves(
            new_curves_temp,self.curve_indices,self.curve_type, self.curve_physical_groups, self.surfaces
        )
        print('Completed removal of duplicate curves.')

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal",1)
        gmsh.model.add(self.name) 
        occ_kernel      = gmsh.model.occ

        # CREATE POINTS
        for i, point in enumerate(self.point_coordinates):
            occ_kernel.addPoint(point[0], point[1], point[2], self.point_mesh_size[i], i + 1)
        print('Created all points.')

        # INSERT BOOLEAN OPERATIONS FOR POINTS

        # CREATE CURVES
        curve_physical_group_indices = []
        for i, curve in enumerate(self.curve_indices):
            if self.curve_type[i] == 0: # LINE
                asdf  =  occ_kernel.addLine(self.curves[curve[0]] + 1, self.curves[curve[1] - 1] + 1, unique_surfaces[i] + 1)
            elif self.curve_type[i] == 1: # CIRCLE ARC
                occ_kernel.addCircleArc(self.curves[curve[0]] + 1, self.curves[curve[1] - 1] + 1, self.curves[curve[1] - 2] + 1, unique_surfaces[i] + 1)
            else:
                pass
            if self.curve_physical_groups[i]:
                curve_physical_group_indices.append(i+1)

        print('Created all curves.')

        # INSERT BOOLEAN OPERATIONS FOR CURVES
        
        # CREATE SURFACES
        surface_physical_group_indices = []
        for i, surface in enumerate(self.surface_indices):
            curve_input = list(self.surfaces[np.arange(surface[0],surface[1])]+1)
            print(curve_input)
            curveloop = occ_kernel.addCurveLoop(curve_input)  # fix the  ccc via Geometry.OCCAutoFix = 0 later
            print(curveloop)
            surface_ind = occ_kernel.addPlaneSurface([curveloop],i+1)

            if self.surface_physical_groups[i]:
                surface_physical_group_indices.append(i+1)
        
        if not surface_physical_group_indices:
            print(surface_physical_group_indices)
            surface_physical_group_indices.append(len(self.surface_indices))
        print('Created all surfaces.')
        print(surface_physical_group_indices)
        # EXECUTE SURFACE BOOLEAN OPERATIONS
        surface_bool_physical_group_indices = []
        for i, parameters in enumerate(self.boolean_parameters):
            print(parameters)
            print([(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_object_indices[i][0],self.boolean_object_indices[i][1])])
            print([(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_tool_indices[i][0],self.boolean_tool_indices[i][1])])
            print(1 + surface_ind + i)
            
            if parameters[0] == 'subtract':
                bool_surf = occ_kernel.cut(
                    [(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_object_indices[i][0],self.boolean_object_indices[i][1])],
                    [(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_tool_indices[i][0],self.boolean_tool_indices[i][1])],
                    # 1 + surface_physical_group_indices[-1] + i + 1,
                    1 + surface_ind + i,
                    removeObject=self.boolean_remove_object[i],
                    removeTool=self.boolean_tool_object[i]
                )
            else:
                raise KeyError('operation has not been integrated yet')
            if self.boolean_surface_physical_groups is not []:
                if self.boolean_surface_physical_groups[i]:
                    # surface_bool_physical_group_indices.append(i + 1 + surface_physical_group_indices[-1])
                    surface_bool_physical_group_indices.append(bool_surf[0][0][1])

        print('Created all boolean surfaces.')

        occ_kernel.synchronize()
 
        # NOTE: Physical groups MUST be added AFTER synchronize()
        # ADD PHYSICAL GROUPS

        curve_counter  = 0
        for i, group in enumerate(self.curve_physical_groups):
            if self.curve_physical_groups[i]:
                curve_counter += 1
                gmsh.model.addPhysicalGroup(1, [curve_physical_group_indices[curve_counter-1]], group[0])
                gmsh.model.setPhysicalName(1, group[0], group[1])  

        surface_counter = 0
        for i, group in enumerate(self.surface_physical_groups):
            if self.surface_physical_groups[i]:
                surface_counter += 1
                gmsh.model.addPhysicalGroup(2, [surface_physical_group_indices[surface_counter-1]], group[0])
                gmsh.model.setPhysicalName(2, group[0], group[1])  

        boolean_surface_counter = 0
        for i, group in enumerate(self.boolean_surface_physical_groups):
            if self.boolean_surface_physical_groups[i]:
                boolean_surface_counter += 1
                gmsh.model.addPhysicalGroup(2, [surface_bool_physical_group_indices[boolean_surface_counter-1]], group[0])
                gmsh.model.setPhysicalName(2, group[0], group[1])  

        # Code block for adding all entities of a dimension to a physical group
        if self.all_points_physical_group:
            pass
        if self.all_curves_physical_group:
            all_curves = occ_kernel.getEntities(dim=1)
            curves = [curve[1] for curve in all_curves]
            gmsh.model.addPhysicalGroup(1, curves, 1000)
            gmsh.model.setPhysicalName(1, 1000, 'All Curves')

        if self.all_surfaces_physical_group:
            pass
        
        # Manual for the initial magnet-disk example for Ru
        if False:
            gmsh.model.addPhysicalGroup(1, list(np.arange(5,20)), 1000)
            gmsh.model.setPhysicalName(1, 1000, 'Interior Curves')

        # should make physical groups all in 1 or 2 arrays, with one being the
        # index, the other being the dimension (can also add names)

        gmsh.model.mesh.generate(2)
        gmsh.write(self.name +  '.msh')

        # NOTE-S:
        # USE getEntities(0) TO EXTRACT POINTS INFO OF MESH
        # USE getNodes(1, ...) TO EXTRACT NODE INFORMATION ON CURVES

        # ------------------------ GETTING POINT INFORMATION ------------------------
        self.gmsh_point_entities = gmsh.model.getEntities(0)
        # we hold this information below to identify what order that gmsh orders its points
        # this can be different from user input b/c of boolean operations & point regeneration
        self.gmsh_order_point_coords = self.reorder_points_to_gmsh(coordinate_system='cartesian')
        self.gmsh_order_point_coords_polar = self.reorder_points_to_gmsh(coordinate_system='polar')

        # ------------------------ GETTING CURVE INFORMATION ------------------------
        self.gmsh_curve_entities = gmsh.model.getEntities(1) # of the form (1, i)
        self.edge_node_coords = []
        self.edge_node_indices = []
        # for edge parametrization
        # nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim=1, tag=5, includeBoundary=True, returnParametricCoord=False)
        for curve in self.gmsh_curve_entities[:]:
            info = gmsh.model.mesh.getNodes(dim=1, tag = curve[1], includeBoundary=True, returnParametricCoord=False)
            self.edge_node_coords.append(np.array(info[1]).reshape((int(len(info[1])/3),3)))
            self.edge_node_indices.append(info[0]) # the node indices along an edge; last two are the start and end nodes

        # nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(dim=1, includeBoundary=True, returnParametricCoord=False)
        # print(nodeTags)
        # print(nodeCoords)

        print('---')
        print(self.edge_node_coords)
        print(self.edge_node_indices)
        print('---')

        self.edge_node_coords = np.array(self.edge_node_coords)

        if coordinate_system is 'polar':
            for i, edge in enumerate(self.edge_node_coords):
                for j, point in enumerate(edge):
                    self.edge_node_coords[i][j] = [
                        np.arctan2(point[1], point[0]),
                        np.linalg.norm(point[:2]),
                        0.0
                    ]
                # sign check for +/- pi to correct parametrization not done here, but in parametrization step

                
        #         if np.pi in self.edge_node_coords[i][-2] or np.pi in self.edge_node_coords[i][-1]:
        #             print(self.edge_node_coords[i])
        #             pi_ind = np.where(abs(np.pi - self.edge_node_coords[i]) < 1e-8) # gives row and column
        #             if np.sign(self.edge_node_coords[i][pi_ind[0], pi_ind[1]]) != np.sign(self.edge_node_coords[i][pi_ind[0]-1, pi_ind[1]]):
        #                 self.edge_node_coords[i][pi_ind[0], pi_ind[1]] -= 2*np.pi
        #                 print('sign mismatch')
        #             print('positive')
        #             print(pi_ind)

        #         if -np.pi in self.edge_node_coords[i][-2] or -np.pi in self.edge_node_coords[i][-1]:
        #             print('negative')
        # print(self.edge_node_coords[-1])
        # exit()


        # if '-nopopup' not in sys.argv:
        if self.popup == True:
            gmsh.fltk.run()
        gmsh.finalize()

    # --------------------- MISCELLANEOUS ---------------------
    def get_coordinates(self, coord_sys='cartesian'):
        num_pts = len(self.point_coordinates)
        
        if coord_sys == 'cartesian':
            self.mesh_nodes = np.array(self.point_coordinates)
        elif coord_sys == 'polar':
            self.mesh_nodes = np.zeros((num_pts, 2))
            for i in range(num_pts):
                self.mesh_nodes[i] = [
                    np.arctan2(self.point_coordinates[i, 1], self.point_coordinates[i, 0]),
                    np.linalg.norm(self.point_coordinates[i, :2]),
                    # 180/np.pi*np.arctan2(self.point_coordinates[i, 1], self.point_coordinates[i, 0])
                ]
        # self.mesh_nodes  = self.mesh_nodes.reshape((num_pts * 2,))
        return self.mesh_nodes

    def reorder_points_to_gmsh(self, coordinate_system='cartesian'): # don't think I need this
        # need to reorder points and edges to align with the final result from gmsh
        # POINTS
        # for points use self.point_coordinates and reorder them
        # print(self.point_coordinates)
        gmsh_order_point_coords = []
        move_origin = 0
        for point in self.gmsh_point_entities:
            point_coord = gmsh.model.mesh.getNodes(
                    dim=point[0], 
                    tag=point[1], 
                    includeBoundary=True, 
                    returnParametricCoord=False
                )[1]

            if np.linalg.norm(point_coord - np.array([0., 0., 0.])) > 1.e-8:
                gmsh_order_point_coords.append(point_coord)
            else:
                move_origin = 1
        
        if move_origin:
            gmsh_order_point_coords.append([0., 0., 0.])

        # print(gmsh_order_point_coords)
        # print(np.array([0., 0., 0.]))
        # origin_ind = np.where(gmsh_order_point_coords[:] == np.array([0., 0., 0.]))
        gmsh_order_point_coords = np.array(gmsh_order_point_coords)

        if coordinate_system is 'polar':
            for i, coords in enumerate(gmsh_order_point_coords):
                gmsh_order_point_coords[i] = [
                    np.arctan2(coords[1], coords[0]),
                    np.linalg.norm(coords[:2]),
                    0.0,
                ]

        return gmsh_order_point_coords
        # we can return the points in the order that GMSH now defines them, in either polar or cartesian
        # need this to relate the point objects to the GMSH defined points

    def test_ffd_edge_parametrization_polar(self, delta, output_type):
        if delta.shape != (int(4 * self.num_ffd_faces), 2):
            raise TypeError('Shape of delta numpy array must be 4 * # FFD faces * 2')
        delta = delta.reshape((2 * delta.shape[0], ))

        # --- FFD TEST
        ordered_coords = self.gmsh_order_point_coords_polar[:,:2].reshape((2 * self.gmsh_order_point_coords_polar.shape[0],))
        point_delta = self.ffd_face_sps_mat.dot(delta)
        new_coords = ordered_coords + point_delta

        # --- EDGE COORDINATE TEST
        new_edge_coords = self.edge_param_sps_mat.dot(new_coords)
        old_edge_coords = self.edge_param_sps_mat.dot(ordered_coords)
        
        if output_type is 'cartesian':
            for i in range(int(len(new_edge_coords)/2)):
                new_edge_coords[2*i:2*i+2] = [
                    new_edge_coords[2*i+1] * np.cos(new_edge_coords[2*i]),
                    new_edge_coords[2*i+1] * np.sin(new_edge_coords[2*i]),
                ]
                old_edge_coords[2*i:2*i+2] = [
                    old_edge_coords[2*i+1] * np.cos(old_edge_coords[2*i]),
                    old_edge_coords[2*i+1] * np.sin(old_edge_coords[2*i]),
                ]

        edge_deltas = new_edge_coords - old_edge_coords
        edge_indices = []
        for i in range(int(len(edge_deltas)/2)):
            edge_indices.extend([i + 1, i + 1])

        return new_edge_coords, old_edge_coords, edge_deltas, np.array(edge_indices) # need to convert to cartesian



# ------------------------------------------- FFD -------------------------------------------
class FFD(object):
    def __init__(self):
        self.entities   = []

        self.vertices       = []

        self.edges          = []
        self.edge_indices   = []
        self.edge_type      = []

        self.faces          = []
        self.face_indices   = []
        self.face_type      = []
        self.ffd_grid_steps = []

        self.blocks          = []
        self.block_indices   = []

    # def add_design_variable(self, dv_name=None, dv1=None, dv1_bounds=None, dv2=None, dv2_bounds=None, dv_type = 'cartesian'):
    #     # dv_type can vary between cartesian and polar
    #     # inputs need to be:
    #     # - design variable directions/axes of interest (so x, y or both) --> corresponds to (r, theta or both) for polar 

    #     pass

    def add_shape_parameter(self, dv_name=None, dv1_bounds=None, dv2_bounds=None, dv_type = 'cartesian'):
        self.dv_dict = {}
        axes = []
        bounds = []

        if dv_type == 'cartesian':
            pass
        elif dv_type == 'polar':
            # dv1_bounds = radius; dv2_bounds = theta
            if dv1_bounds is not None:
                axes.append('r')
                bounds.append(dv1_bounds)

            if dv2_bounds is not None:
                axes.append('theta')
                bounds.append(dv2_bounds)

        self.dv_dict['name'] = dv_name
        self.dv_dict['axes'] = axes
        self.dv_dict['bounds'] = bounds
        print(self.dv_dict)

        return(self.dv_dict)


    def link_design_variable(self, face, design_variable):
        pass

    def add_ffd_entity(self, entity, u_p=10, v_p=10):
        self.entities.append(entity)
        self.ffd_grid_steps.append([u_p,v_p])

    def create_ffd_grid(self, ffd_grid_steps=None):
        if ffd_grid_steps is None:
            ffd_grid_steps = []

        for i, step in enumerate(ffd_grid_steps):

            edges           = self.faces[np.arange(self.face_indices[i,0],self.face_indices[i,1])]
            face_vertices_ind   = []

            for j, edge in enumerate(edges):
                vertices    = self.edges[np.arange(self.edge_indices[j,0],self.edge_indices[j,1])]
                face_vertices_ind.extend(vertices)

            face_vertices_ind_set   = list(set(face_vertices_ind))
            print(face_vertices_ind)

            face_vertices   = [list(self.vertices[i]) for i in face_vertices_ind]
            print(face_vertices)


    def define_vertex(self, x, y, z):
        self.vertices.append([x, y, z])
        return len(self.vertices) - 1

    def define_edge(self, children, edge_type):
        self.edges.extend([vertex.id for vertex in children])
        if self.edge_indices == []:
            self.edge_indices.append([0, len(children)])
        else:
            self.edge_indices.append([
                self.edge_indices[-1][-1],
                self.edge_indices[-1][-1]  + len(children),
            ])
        return len(self.edge_indices) - 1

    def define_face(self, children):
        self.faces.extend([edge.id for edge in children])
        if self.face_indices == []:
            self.face_indices.append([0, len(children)])
        else:
            self.face_indices.append([
                self.face_indices[-1][-1],
                self.face_indices[-1][-1] + len(children)
            ])

        return len(self.face_indices) - 1

    def define_block(self):
        pass

    def assemble(self):
        for entity in self.entities:
            print(1)
            entity.assemble(self)

        self.vertices       = np.array(self.vertices)

        self.edges          = np.array(self.edges)
        self.edge_indices   = np.array(self.edge_indices)
        self.edge_type      = np.array(self.edge_type)

        self.faces          = np.array(self.faces)
        self.face_indices   = np.array(self.face_indices)
        self.face_type      = np.array(self.face_type)
        self.ffd_grid_steps = np.array(self.ffd_grid_steps)

        # self.blocks         = np.array(self.blocks)
        # self.block_indices  = np.array(self.block_indices)
        # self.block_type     = np.array(self.block_type)

        # removing duplicate vertices
        self.vertices, self.edges    = remove_duplicate_points(
            points=self.vertices, 
            point_mesh_size=None, 
            curves=self.edges,
        )

        # removing duplicate edges
        self.edges, self.edge_indices, self.faces, unique_faces = remove_duplicate_curves(
            curves=self.edges, 
            curve_indices=self.edge_indices, 
            curve_type = None, 
            surfaces=self.faces,
        )

        # print('---')
        # print('vertices:', self.vertices)
        # print('edges:', self.edges)
        # print('edge indices:', self.edge_indices)
        # print('faces:', self.faces)
        # print('face indices:', self.face_indices)
        # print('---')

        # ffd_grid    = self.create_ffd_grid(self.ffd_grid_steps)


    # NOTE: remove the FFD class
    # - duplicates are ok, the assemble step is not necessary  (like it is for Mesh)