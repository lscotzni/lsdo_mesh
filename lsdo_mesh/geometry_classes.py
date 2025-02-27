from re import L
import gmsh
import numpy as np
import scipy as sp
from scipy.linalg import norm
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix
import sys
import copy

import csdl
from csdl import Model

from lsdo_mesh.mesh_ffd_objects import *
from lsdo_mesh.remove_duplicates import remove_duplicate_points, remove_duplicate_curves
from lsdo_mesh.csdl_mesh_models import EdgeUpdateModel, ShapeParameterModel

from lsdo_mesh.geometry_operations import rotate
import os

# ------------------------------------------- MESH -------------------------------------------
class Mesh(object):

    def __init__(self, name='Mesh', dim=2, popup=False, rotation_angles=[0]):

        self.name                   = name
        self.popup                  = popup
        self.rotation_angles        = rotation_angles
        self.dim                    = dim
        
        self.top_entities           = []

        self.point_coordinates      = []
        self.point_mesh_size        = []
        self.point_rotate_instance  = []

        self.curves                 = []
        self.curve_indices          = []
        self.curve_type             = []

        self.surfaces               = []
        self.surface_indices        = []
        self.surface_curve_loops    = []
        self.surface_type           = []
        self.surface_colors_all     = []

        self.boolean_operations     = []
        self.boolean_entities       = []
        self.boolean_object_indices = []
        self.boolean_tool_indices   = []
        self.boolean_remove_object  = []
        self.boolean_tool_object    = []
        self.boolean_parameters     = []
        self.boolean_surface_colors_all = []

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

        self.mesh_points_instances  = [] # holds points of all instances of rotations
        # length of above variable is equal to number of rotation instances
        self.ffd_cp_instances       = []
        self.edge_nodes_instances   = []
        self.gmsh_order_point_coords_instances  = []

        # variables below are for point indices stored 
        self.point_node_indices_to_file             = False
        self.points_for_node_extraction             = []

        self.point_node_coords_to_file              = False
        self.points_for_coord_extraction            = []
        
        self.extracted_node_instances               = []
        self.extracted_node_coord_instances         = []

    # --------------------- GEOMETRY/MESH ---------------------
    def add_entity(self, entity=None):
        # coordinates is the indicator for either polar or cartesian for full entity
        self.top_entities.append(entity)

    def add_point(self, x, y, z, mesh_size, rotate_instance):
        self.point_coordinates.append([x, y, z])
        self.point_mesh_size.append(mesh_size)
        self.point_rotate_instance.append(rotate_instance)
        return len(self.point_coordinates) - 1

    def add_curve(self, children, curve_type, physical_group=False):
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

        if type(physical_group) is not tuple:
            if physical_group is not False:
                raise KeyError('Physical groups must be tuples denoted as (int, string).')
        self.curve_physical_groups.append(physical_group)
        
        return len(self.curve_type) - 1

    def add_surface(self, children, curve_loop_lengths, physical_group=False, color=False):
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
                raise TypeError('Physical groups must be tuples denoted as (int, string).')
        self.surface_physical_groups.append(physical_group)
        self.surface_curve_loops.append(curve_loop_lengths)
        self.surface_colors_all.append(color)

        return len(self.surface_indices) - 1

    def add_boolean(self, children, properties, physical_group=False, color=False):
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
        self.boolean_surface_colors_all.append(color)

        return len(self.boolean_operations) - 1

    def add_physical_group(self, entities):

        if not all([isinstance(entities[0], entity) for entity in entities[1:]]):
            raise KeyError('Inputs must be the same type.')

    def add_all_entities_to_physical_group(self, geometry_type=None):

        if geometry_type == 'points':
            raise KeyError('Not implemented with points yet.')
        elif geometry_type == 'curves':
            dim = 1
            self.all_curves_physical_group = True
        elif geometry_type == 'surfaces':
            raise KeyError('Not implemented with surfaces yet.')

    # --------------------- FFD ---------------------
    def add_face(self, face):
        self.ffd_faces.append(face)

    # --------------------- ASSEMBLE ---------------------
    def assemble(self, coordinate_system='cartesian'):
        self.assemble_mesh(coordinate_system=coordinate_system)

        self.get_coordinates(coord_sys=coordinate_system) # returns self.mesh_nodes for entire mesh

        # PARAMETRIZATION STEPS
        self.assemble_shape_parameter_parametrization(coordinate_system=coordinate_system)

        self.assemble_ffd_control_points(coordinate_system=coordinate_system)

        self.assemble_ffd_parametrization(coordinate_system=coordinate_system)

        self.assemble_edge_parametrization(coordinate_system=coordinate_system)

        self.assemble_mesh_parametrization(coordinate_system=coordinate_system) # ---- Ru's work here ----
        # long term, line above contains Ru's FEniCS mesh deformation/regeneration model
        # for short term, we will separate it

        self.return_ffd_parameters(coordinate_system=coordinate_system)

        self.create_csdl_model() # this contains the ffd face/edge & mesh movement csdl model classes
        # spits out the csdl variable containing mesh coordinates
        print('Completed lsdo_mesh assemble method.')
        return self.mesh_models, self.shape_parameter_model

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
        self.point_rotate_instance      = np.array(self.point_rotate_instance)
        self.point_physical_groups      = np.array(self.point_physical_groups)

        self.curves                     = np.array(self.curves)
        self.curve_indices              = np.array(self.curve_indices)
        self.curve_type                 = np.array(self.curve_type)
        # self.curve_physical_groups      = np.array(self.curve_physical_groups)

        self.surfaces                   = np.array(self.surfaces)
        self.surface_indices            = np.array(self.surface_indices)
        # self.surface_curve_loops        = np.array(self.surface_curve_loops)
        # not turning self.surface_curve_loops into np.array b/c it contains lists of different lengths
        self.surface_type               = np.array(self.surface_type)
        # self.surface_physical_groups    = np.array(self.surface_physical_groups)

        self.boolean_operations         = np.array(self.boolean_operations)
        self.boolean_entities           = np.array(self.boolean_entities)
        self.boolean_object_indices     = np.array(self.boolean_object_indices)
        self.boolean_tool_indices       = np.array(self.boolean_tool_indices)
        self.boolean_remove_object      = np.array(self.boolean_remove_object)
        self.boolean_tool_object        = np.array(self.boolean_tool_object)
        self.boolean_parameters         = np.array(self.boolean_parameters)

        # removing duplicate points
        print('Starting removal of duplicate points.')
        self.point_coordinates, self.point_mesh_size, self.point_rotate_instance, new_curves_temp = remove_duplicate_points(
            self.point_coordinates, self.point_mesh_size, self.point_rotate_instance, self.curves, 
            )
        print('Completed removal of duplicate points.')

        self.num_points     = self.point_coordinates.shape[0]

        # removing duplicate curves
        # print('Starting removal of duplicate curves.')
        self.curves, self.curve_indices, self.curve_type, self.curve_physical_groups, self.surfaces, unique_surfaces = remove_duplicate_curves(
            new_curves_temp,self.curve_indices,self.curve_type, self.curve_physical_groups, self.surfaces
        )
        # print('Completed removal of duplicate curves.')
        aaa = []
        bbb = []
        for a, angle in enumerate(self.rotation_angles):
            self.instance   = a
            point_rotate_instances = list(self.point_rotate_instance)

            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal",1)
            gmsh.model.add(self.name) 
            occ_kernel      = gmsh.model.occ

            # CREATE POINTS
            for i, point in enumerate(self.point_coordinates):
                occ_kernel.addPoint(point[0], point[1], point[2], self.point_mesh_size[i], i + 1)
                if point_rotate_instances[i]:
                    occ_kernel.rotate([(0, i+1)], 0, 0, 0, 0, 0, 1, angle)
            print('Created all points.')

            # INSERT BOOLEAN OPERATIONS FOR POINTS

            # CREATE CURVES
            curve_physical_group_indices = []
            for i, curve in enumerate(self.curve_indices):
                if self.curve_type[i] == 0: # LINE
                    asdf  =  occ_kernel.addLine(self.curves[curve[0]] + 1, self.curves[curve[1] - 1] + 1, unique_surfaces[i] + 1)
                elif self.curve_type[i] == 1: # CIRCLE ARC
                    occ_kernel.addCircleArc(self.curves[curve[0]] + 1, self.curves[curve[1] - 1] + 1, self.curves[curve[1] - 2] + 1, unique_surfaces[i] + 1)
                # else:
                #     pass
                if self.curve_physical_groups[i]:
                    curve_physical_group_indices.append(i+1)

            print('Created all curves.')

            # INSERT BOOLEAN OPERATIONS FOR CURVES

            # CREATE CURVE LOOPS
            
            # CREATE SURFACES
            surface_physical_group_indices = []
            surface_color_indices   = []
            surface_colors          = []
            for i, surface in enumerate(self.surface_indices):
                # for ... in <indicator for number of curve loops here>
                curve_loop_lengths = self.surface_curve_loops[i]
                # print(curve_loop_lengths)
                
                gmsh_curve_loops = []
                loop_counter = 0
                for j, loop_size in enumerate(curve_loop_lengths):
                    num_curves_in_loop = loop_size
                    curve_input = list(self.surfaces[np.arange(surface[0] + loop_counter, surface[0] + loop_counter + loop_size)]+1)
                    # print(curve_input)
                    curveloop = occ_kernel.addCurveLoop(curve_input)  # fix the  ccc via Geometry.OCCAutoFix = 0 later
                    gmsh_curve_loops.append(curveloop)
                    loop_counter += loop_size
                # surface_ind = occ_kernel.addPlaneSurface([curveloop],i+1)
                surface_ind = occ_kernel.addPlaneSurface(gmsh_curve_loops,i+1)
                # print(gmsh_curve_loops)
                # exit()
                if self.surface_physical_groups[i]:
                    # surface_physical_group_indices.append(i+1)
                    surface_physical_group_indices.append(surface_ind)

                if self.surface_colors_all[i]:
                    surface_color_indices.append(surface_ind)
                    surface_colors.append(self.surface_colors_all[i])
            
            if not surface_physical_group_indices:
                surface_physical_group_indices.append(len(self.surface_indices))
            print('Created all surfaces.')

            # EXECUTE SURFACE BOOLEAN OPERATIONS
            bool_surface_physical_group_indices = []
            bool_surface_color_indices  = []
            bool_surface_colors         = []
            for i, parameters in enumerate(self.boolean_parameters):
                
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
                        # bool_surface_physical_group_indices.append(i + 1 + surface_physical_group_indices[-1])
                        bool_surface_physical_group_indices.append(bool_surf[0][0][1])
                
                if self.boolean_surface_colors_all[i]:
                    bool_surface_color_indices.append(bool_surf[0][0][1])
                    bool_surface_colors.append(self.boolean_surface_colors_all[i])

            print('Created all boolean surfaces.')

            occ_kernel.synchronize()
    
            # NOTE: Physical groups and color setting MUST be added AFTER synchronize()

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
                    gmsh.model.addPhysicalGroup(2, [bool_surface_physical_group_indices[boolean_surface_counter-1]], group[0])
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

            # SET COLORS
            # SURFACE COLORS
            for i in range(len(surface_colors)):
                gmsh.model.setColor(
                    [(2,surface_color_indices[i])],
                    surface_colors[i][0],
                    surface_colors[i][1],
                    surface_colors[i][2],
                )

            # BOOLEAN SURFACE COLORS
            for i in range(len(bool_surface_colors)):
                gmsh.model.setColor(
                    [(2,bool_surface_color_indices[i])],
                    bool_surface_colors[i][0],
                    bool_surface_colors[i][1],
                    bool_surface_colors[i][2],
                )

            gmsh.model.mesh.generate(2)
            gmsh.write(self.name + '_{}.msh'.format(str(a+1)))

            # NOTE-S:
            # USE getEntities(0) TO EXTRACT POINTS INFO OF MESH
            # USE getNodes(1, ...) TO EXTRACT NODE INFORMATION ON CURVES

            # ------------------------ GETTING POINT INFORMATION ------------------------
            # Getting point entities and appending the point instances to a global array
            self.gmsh_point_entities    = gmsh.model.getEntities(0) # USER-DEFINED POINTS of the form (0, i)
            self.mesh_points_instances.append(np.array(
                self.reorder_points_to_gmsh(coordinate_system=coordinate_system)[0][:,:2].reshape((self.num_points * 2, ))
            ))
            self.gmsh_order_point_coords_instances.append(np.array(
                self.reorder_points_to_gmsh(coordinate_system='cartesian')[0]
            ))
            # We need edge and point information for each mesh because of how parametrization/arctan2 works around pi
            # 
            if self.points_for_node_extraction: # ONLY WANT IF USER CALLS self.get_node_indices & ADDS POINTS
                self.gmsh_point_nodes           = gmsh.model.mesh.getNodes(dim=0) # node information for POINTS; only needed in this case
                extracted_point_node_indices    = self.extract_point_node_indices(self.points_for_node_extraction)
                self.extracted_node_instances.append(extracted_point_node_indices)
            
            if self.points_for_coord_extraction:
                self.gmsh_point_nodes           = gmsh.model.mesh.getNodes(dim=0) # node information for POINTS; only needed in this case
                extracted_point_node_coords     = self.extract_point_node_coords(self.points_for_coord_extraction)
                self.extracted_node_coord_instances.append(extracted_point_node_coords) 
            
            # ------------------------ GETTING CURVE INFORMATION ------------------------
            self.gmsh_curve_entities = gmsh.model.getEntities(1) # of the form (1, i)
            edge_node_coords = []
            # for edge parametrization
            for curve in self.gmsh_curve_entities[:]:
                info = gmsh.model.mesh.getNodes(dim=1, tag = curve[1], includeBoundary=True, returnParametricCoord=False)
                edge_node_coords.append(np.array(info[1]).reshape((int(len(info[1])/3),3)))
            
            edge_node_coords = np.array(edge_node_coords) # specify dtype to avoid warning
            
            if coordinate_system == 'polar':
                for i, edge in enumerate(edge_node_coords):
                    for j, point in enumerate(edge):
                        edge_node_coords[i][j] = [
                            np.arctan2(point[1], point[0]),
                            np.linalg.norm(point[:2]),
                            0.0
                        ]

            # print(edge_node_coords)
            # print(edge_node_coords[0])
            # print(edge_node_coords[1])
            # exit()
            
            self.edge_nodes_instances.append(edge_node_coords)

            # STUFF THAT DOESN'T CHANGE WITH MESH INSTANCES (BE CAREFUL, THIS MAY BE A CONCERING ASSUMPTION)
            if a == 0:
                self.gmsh_order_point_node_ind = self.reorder_points_to_gmsh(coordinate_system='cartesian')[1] # NEED ONLY ONCE
                self.edge_node_indices = []
                for curve in self.gmsh_curve_entities[:]:
                    info = gmsh.model.mesh.getNodes(dim=1, tag = curve[1], includeBoundary=True, returnParametricCoord=False)
                    self.edge_node_indices.append(info[0]) # the node indices along an edge; last two are the start and end nodes
            print('---')
            # print(gmsh.model.mesh.getNodes(dim=1, includeBoundary=True, returnParametricCoord=False)[1])
            aaa.append(gmsh.model.mesh.getNodes(dim=1, includeBoundary=True, returnParametricCoord=False)[0])
            bbb.append(gmsh.model.mesh.getNodes(dim=1, includeBoundary=True, returnParametricCoord=False)[1])
            # print(len(aaa[0]))
            print('---')

            # if '-nopopup' not in sys.argv:
            if self.popup == True:
                gmsh.fltk.run()

            gmsh.finalize()

            os.system('python3 msh2xdmf.py -d 2 ' + self.name + '_{}.msh'.format(str(a+1)))

        # print(aaa[1] - aaa[0])
        # coord_diff = bbb[1] - bbb[0]
        # asdf = np.where(coord_diff != 0.)
        # print(coord_diff[asdf])
        # print(asdf)
        # self.mesh_points_instances = np.array(self.mesh_points_instances)

    def assemble_shape_parameter_parametrization(self, coordinate_system='cartesian'):
        sparse_val, sparse_row, sparse_col = [], [], []
        column_counter = 0

        self.shape_parameter_list    = []
        self.axis_list               = []
        self.parameter_type_list     = [] # holds strings (constant & linear)
        self.parameter_index_list    = [] # holds 1 for constant and 2 for linear

        for i, face in enumerate(self.ffd_faces):

            # sparse matrix formatting:
            # for x or theta, we want the odd rows for each set of 8; in python, this is the even rows (0,2,4,6)
            # for y or r, we want the even rows for each set of 8; in python, this is the odd rows (1,3,5,7)
            
            all_face_shape_parameters = vars(face)['parameters']
            for j, parameters in enumerate(all_face_shape_parameters):
                # Adding parameters to lists in order of definition
                self.shape_parameter_list.append(parameters[0] + '_sp')
                self.axis_list.append(parameters[1])
                self.parameter_type_list.append(parameters[2])

                # format for parameters: [name, axis, def_type]
                if parameters[1] in ['x', 'theta']:
                    dim = 0 # applied to first coordinate
                    sparse_row.extend(list(np.arange(8 * i, 8 * (i+1) - 1, 2)))
                elif parameters[1] in ['y', 'r']:
                    dim = 1 # applied to second coordinate
                    sparse_row.extend(list(np.arange(8 * i + 1, 8 * (i+1), 2)))
                
                if parameters[2] == 'constant':
                    self.parameter_index_list.append(1) # adding 1 for constant

                    col = 1
                    u   = 1
                    # contains 1 entry in vector
                    sparse_val.extend([1] * 4)
                    sparse_col.extend([column_counter] * 4)

                elif parameters[2] == 'linear':
                    self.parameter_index_list.append(2) # adding 2 for linear

                    col = 2
                    u   = 1/2  # use negative for one way, positive for the other
                    # contains 2 entries in vector
                    if dim == 0:
                        sparse_col.extend([column_counter, column_counter + 1] * 2)
                        vec = [-1/2, 1/2, -1/2, 1/2]
                        #  positive direction for x and theta differ, so sign change is made here
                        if parameters[1] == 'x':
                            sparse_val.extend(vec)
                        elif parameters[1] == 'theta':
                            sparse_val.extend([1 * entry for entry in vec])
                    elif dim == 1:
                        sparse_val.extend([-1/2, -1/2, 1/2, 1/2])
                        sparse_col.extend([
                            column_counter, column_counter, column_counter + 1, column_counter + 1,
                        ])

                column_counter += col

        # sparse matrix produced here is of shape:
        # (num FFD faces * 8, num shape parameters)

        self.shape_param_sps_mat = csc_matrix(
            (sparse_val, (sparse_row, sparse_col)),
            shape=((len(self.ffd_faces) * 8, max(sparse_col) + 1))
        )

    def assemble_ffd_control_points(self, coordinate_system='cartesian'):
        self.num_ffd_faces          = len(self.ffd_faces)
        self.num_ffd_face_coords    = 4 * self.dim

        self.ffd_face_control_pts   = np.zeros(
            (self.num_ffd_face_coords * self.num_ffd_faces, )
        )

        # LOOP GETTING FFD FACE CP COORDINATES FOR FIRST GEOMETRY (NOT ROTATED)
        for face_ind, face in enumerate(self.ffd_faces):
            face_coords     = np.array(face.return_coordinates(coordinate_system)) # reshape((4 * dim)) to vectorize

            # physical coordinates of ffd face vertices
            # may be good to hold this information in the vertex object
            # i.e. x1 = [x1_min, x1_max]; x2 = [x2_min, x2_max]
            P00 = [np.min(face_coords[:,0]), np.min(face_coords[:,1])]
            P01 = [np.min(face_coords[:,0]), np.max(face_coords[:,1])]
            P10 = [np.max(face_coords[:,0]), np.min(face_coords[:,1])]
            P11 = [np.max(face_coords[:,0]), np.max(face_coords[:,1])]

            face_coords = np.array([P00, P10, P01, P11])

            start       = self.num_ffd_face_coords * face_ind
            end         = start + self.num_ffd_face_coords
            # in the order of P00, P10, P01, P11
            self.ffd_face_control_pts[start:end] = face_coords.reshape((self.dim * len(face_coords), ))

        # self.ffd_cp_instances.append(self.ffd_face_control_pts)
        self.ffd_cp_instances.append(self.ffd_face_control_pts)

        '''
        WE ARE CREATING COORDINATE INSTANCES OF ALL OF THE FFD CONTROL POINTS FOR ROTATIONS
        WE NEED TO CHECK WHICH ONES GET ROTATED SO IN THE SECOND LOOP
        FOR POINTS THAT LIE WITHIN FFD FACES THAT GET ROTATED, WE MANUALLY ROTATE THOSE POINTS
        '''

        for a, angle in enumerate(self.rotation_angles):
            if a == 0:
                continue

            ffd_face_control_pts = []
            ffd_face_control_pts.extend(self.ffd_face_control_pts)
            # ffd_face_control_pts = self.ffd_face_control_pts[:]
            for face_ind, face in enumerate(self.ffd_faces):
                embedded_points = vars(face)['embedded_points']
                face_dummy = copy.copy(face)

                if all([vars(embedded_points[i])['rotate_instance'] for i in range(len(embedded_points))]):
                    face_copy = copy.copy(face_dummy)
                    control_points  = vars(face_copy)['face_vertices']
                    # print('before:', control_points)
                    control_points  = rotate(
                        entities=control_points,
                        angle=[0., 0., angle],
                    )
                    vars(face_copy)['face_vertices'] = control_points
                    # print('after:', control_points)
                
                    face_coords     = np.array(face_copy.return_coordinates(coordinate_system)) # reshape((4 * dim)) to vectorize
                    # face_coords     = [control_points[i].return_coordinates(coordinate_system) for i in range(4)]
                    # face_coords     = np.array(face_coords)

                    # physical coordinates of ffd face vertices
                    # may be good to hold this information in the vertex object
                    # i.e. x1 = [x1_min, x1_max]; x2 = [x2_min, x2_max]
                    P00 = [np.min(face_coords[:,0]), np.min(face_coords[:,1])]
                    P01 = [np.min(face_coords[:,0]), np.max(face_coords[:,1])]
                    P10 = [np.max(face_coords[:,0]), np.min(face_coords[:,1])]
                    P11 = [np.max(face_coords[:,0]), np.max(face_coords[:,1])]

                    face_coords = np.array([P00, P10, P01, P11])

                    start       = self.num_ffd_face_coords * face_ind
                    end         = start + self.num_ffd_face_coords
                    # in the order of P00, P10, P01, P11
                    ffd_face_control_pts[start:end] = face_coords.reshape((self.dim * len(face_coords), ))
                elif not all([vars(embedded_points[i])['rotate_instance'] for i in range(len(embedded_points))]):
                    if any([vars(embedded_points[i])['rotate_instance'] for i in range(len(embedded_points))]):
                        raise KeyError('All embedded points in FFD face must rotate together.')

            # self.ffd_cp_instances.append(np.array(ffd_face_control_pts))
            self.ffd_cp_instances.append(ffd_face_control_pts)
        self.ffd_cp_instances = np.array(self.ffd_cp_instances)
        # exit()

        # # CHECK FFD CP AND DIFFERENCE FOR EACH ROTATION INSTANCE
        # for i in range(len(self.rotation_angles)):
        #     print(self.ffd_cp_instances[i])

        # for i in range(len(self.rotation_angles) - 1):
        #     print(self.ffd_cp_instances[i+1] - self.ffd_cp_instances[i])

        # exit()
        
    def assemble_ffd_parametrization(self, coordinate_system='cartesian'):
        
        print(' ============ ASSEMBLING FFD FACE PARAMETRIZATION ============ ')

        self.num_ffd_faces   = len(self.ffd_faces)
        num_ffd_face_coords = 4 * self.dim # number of coordinate components stored (2D is x1, x2, etc.)'

        # num col = # of FFD FACES * 4 * 2 
        # P_new = SPS_MAT * (V - V0) + P0

        self.ffd_face_sps_mat_list = []
        # START ROTATION HERE
        for a, angle in enumerate(self.rotation_angles):
            sparse_row, sparse_col, sparse_val = [], [], []
            print('angle instance:', a)

            ffd_face_control_pts = self.ffd_cp_instances[a]

            for face_ind, face in enumerate(self.ffd_faces):
                start   = self.num_ffd_face_coords * face_ind
                end     = start + self.num_ffd_face_coords
                
                # need to extract P00 and P11 from self.ffd_face_control_pts
                P00     = [ffd_face_control_pts[start], ffd_face_control_pts[start + 1]]
                P11     = [ffd_face_control_pts[end - 2], ffd_face_control_pts[end - 1]]

                # loop for the number of children within the face
                embedded_points = vars(face)['embedded_points'] # need to apply rotation here for points with rotated instances

                # CHECKING THE DIMENSION OF THE ENTITY IN THE FFD FACE FOR PARAMETRIZATION
                dim_0_diff_norm = abs(P11[0] - P00[0])
                dim_1_diff_norm = abs(P11[1] - P00[1])
                if len(embedded_points) == 1: # THIS IS A SINGLE POINT WITHIN THE FFD FACE
                    embedded_dim = 0
                # elif (P00[0] == P11[0]) and (P00[1] != P11[1]): # THIS IS A CURVE CONSTANT IN DIMENSION 1 (x OR theta)
                #     embedded_dim = 1
                #     const_dim = 0
                # elif (P00[0] != P11[0]) and (P00[1] == P11[1]): # THIS IS A CURVE CONSTANT IN DIM 2 (y OR r)
                #     embedded_dim = 1
                #     const_dim = 1
                elif (dim_0_diff_norm < 1e-6) and (dim_1_diff_norm > 1e-6): # THIS IS A CURVE CONSTANT IN DIMENSION 1 (x OR theta)
                    embedded_dim = 1
                    const_dim = 0
                elif (dim_0_diff_norm > 1e-6) and (dim_1_diff_norm < 1e-6): # THIS IS A CURVE CONSTANT IN DIM 2 (y OR r)
                    embedded_dim = 1
                    const_dim = 1
                else:
                    embedded_dim = 2

                # adjust face coordinates for discontinuity at "-" x-axis
                theta_discontinuity = 0
                if P00[0] < 0 and P11[0] > 0 and abs(P11[0] - P00[0]) > np.pi:
                    theta_discontinuity = 1
                    # P00[0], P11[0] = P11[0], P00[0] + 2*np.pi # ORDER SWAPPED HERE BECAUSE OF 2pi ADDITION
                    P00[0], P11[0] = P11[0], P00[0] + 2*np.pi # ORDER SWAPPED HERE BECAUSE OF 2pi ADDITION

                for point_ind, point in enumerate(embedded_points):
                    if a != 0:
                        if all([vars(embedded_points[i])['rotate_instance'] for i in range(len(embedded_points))]):
                            point = rotate(
                                entities=[point],
                                angle=[0., 0., angle]
                            )[0]
                    point_coords_to_reorder = np.array(point.return_coordinates(output_type='cartesian'))
                    point_coords            = point.return_coordinates(coordinate_system)[:2]
        
                    # might be good to keep this comparison in cartesian coordinates
                    # polar coordinates are tough b/c np.arctan2() returns angle in range [-pi, pi]
                    # pi and -pi represent the same thing in our case but the sign will break the comparison
                    # print(embedded_points)
                    # print(point_coords_to_reorder)
                    # print(self.gmsh_order_point_coords_instances[a])
                    index = np.where(
                        np.linalg.norm(
                            point_coords_to_reorder - self.gmsh_order_point_coords_instances[a],
                            axis=1
                        ) < 1.e-6
                    )
                    # print('index:', index)

                    # ADD PI-CONDITION 
                    if theta_discontinuity == 1:
                        if point_coords[0] < 0:
                            point_coords[0] += 2*np.pi
                            
                    if embedded_dim == 0:
                        [u, v] = [1, 1]
                    elif embedded_dim == 1:
                        if const_dim == 0:
                            [u,v] = [1, (point_coords[1] - P00[1]) / (P11[1] - P00[1])]
                        elif const_dim == 1:
                            [u,v] = [(point_coords[0] - P00[0]) / (P11[0] - P00[0]), 1]
                    else:
                        [u, v] = [
                            (point_coords[0] - P00[0]) / (P11[0] - P00[0]),
                            (point_coords[1] - P00[1]) / (P11[1] - P00[1]),
                        ]

                    # if np.any([i < 0 for i in [u,v]])  or np.any([i > 1 for i in [u,v]]):
                    #     print('outside of bounds of 0 & 1')

                    for i in range(4):
                        sparse_row.extend(
                            np.arange(self.dim*index[0][0], self.dim*index[0][0] + self.dim)
                        )
                    sparse_col.extend(np.arange(start, end, dtype=int))
                    sparse_val.extend([
                        (1 - u) * (1 - v),  # P00
                        (1 - u) * (1 - v),  # P00
                        u * (1 - v),        # P10
                        u * (1 - v),        # P10
                        (1 - u) * v,        # P01
                        (1 - u) * v,        # P01
                        u * v,              # P11
                        u * v,              # P11
                    ])
                
            # exit()
            if a == 0:
                self.ffd_face_sps_mat = csc_matrix(
                    (sparse_val, (sparse_row, sparse_col)),
                    shape=(self.num_points * self.dim, self.num_ffd_face_coords * self.num_ffd_faces)
                )
            
            ffd_face_sps_mat = csc_matrix(
                (sparse_val, (sparse_row, sparse_col)),
                shape=(self.num_points * self.dim, self.num_ffd_face_coords * self.num_ffd_faces)
            )
            self.ffd_face_sps_mat_list.append(ffd_face_sps_mat)

    def assemble_edge_parametrization(self, coordinate_system='cartesian'):
        # print(' ============ ASSEMBLING EDGE PARAMETRIZATION ============ ')

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

        num_edge_nodes = max([max(node) for node in self.edge_node_indices])

        [num_points, dim] = self.mesh_nodes.shape
        num_edges = len((self.edge_node_indices))

        self.edge_param_sps_mat_list = []
        # EACH MESH INSTANCE HAS NEW EDGE PARAMETRIZATION
        for a, angle in enumerate(self.rotation_angles): 
            sparse_row, sparse_col, sparse_val = [], [], []

            edge_node_coords = self.edge_nodes_instances[a]
            
            # IDENTITY MAP FOR NODES CORRESPONDING TO USER-DEFINED MESH POINTS
            # KEY NOTE HERE: USING ORDER OF GMSH POINTS
            for i, node_ind in enumerate(self.gmsh_order_point_node_ind):
                sparse_row.extend([self.dim * i, self.dim * i + 1])
                sparse_col.extend([self.dim * i, self.dim * i + 1])
                sparse_val.extend([1.0, 1.0])
            
            # use self.edge_boundary_nodes and self.edge_node_indices
            for i, edge in enumerate(self.edge_node_indices):
                # edge_nodes = edge_node_coords[i]
                edge_nodes = copy.deepcopy(edge_node_coords[i])
                # exit()
                first_edge = edge_nodes[0]
                # print(edge_nodes)
                # exit()

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

                start = np.where(self.gmsh_order_point_node_ind == int(edge[-2]))[0][0]
                end = np.where(self.gmsh_order_point_node_ind == int(edge[-1]))[0][0]
                
                num_internal_nodes = len(edge) - 2

                # WORK-AROUND FOR THE DISCONTINUITY AT THETA=PI
                theta_discontinuity = 0
                min_angle = min(edge_nodes[-2][0], edge_nodes[-1][0])
                max_angle = max(edge_nodes[-2][0], edge_nodes[-1][0])
                if abs(max_angle - min_angle) > np.pi and np.sign(max_angle) != np.sign(min_angle):
                    theta_discontinuity = 1
                    if edge_nodes[-2,0] == min_angle:
                        edge_nodes[-2,0] += 2*np.pi
                    elif edge_nodes[-1,0] == min_angle:
                        edge_nodes[-1,0] += 2*np.pi
                    # min_angle, max_angle = max_angle, min_angle + 2*np.pi
                    # edge_nodes[-2, 0] = min_angle
                    # edge_nodes[-1, 0] = max_angle

                u = []
                for j in range(num_internal_nodes):
                    if theta_discontinuity == 1 and var_dim == 0:
                        if edge_nodes[j, var_dim] < 0:
                            edge_nodes[j, var_dim] += 2*np.pi
                    u.append(
                        (edge_nodes[j, var_dim] - edge_nodes[-2, var_dim]) / \
                        (edge_nodes[-1, var_dim] - edge_nodes[-2, var_dim])
                    )

                    sparse_row.extend(
                        [int(2 * (edge[j] - 1))] * 2 + [int(2 * (edge[j] - 1) + 1)] * 2,
                    )
                    sparse_col.extend([
                        2 * start,
                        2 * end,
                        2 * (start) + 1, 
                        2 * (end) + 1,
                    ])
                    sparse_val.extend([
                        (1 - u[j]),
                        u[j], 
                        (1 - u[j]),
                        u[j],
                    ])
            
            edge_param_sps_mat = csc_matrix(
                (sparse_val, (sparse_row, sparse_col)),
                shape=(int((num_edge_nodes + 1) * self.dim), int((num_points) * self.dim)),
            )
            self.edge_param_sps_mat_list.append(edge_param_sps_mat)
   
    def assemble_mesh_parametrization(self, coordinate_system='cartesian'):
        pass

    def create_csdl_model(self):

        # STEPS:
        # 1. OPTIMIZER OUTPUTS UPDATED SHAPE PARAMETER VALUES
        # 2. APPLY MAT-VEC PRODUCT BETWEEN SHAPE PARAMETER SPARSE MAT AND SHAPE PARAMETER DELTAS;
        #       - THEN ADD TO ORIGINAL FFD CPs (FOUND IN ffd_cp_instances)
        # 3. APPLY MAT-VEC PRODUCT BETWEEN FFD SPARSE MAT AND NEW FFD CPs
        #       - ADD TO ORIGINAL POINT COORDINATES (FOUND IN mesh_points_instances)
        # 4. APPLY MAT-VEC PRODUCT BETWEEN EDGE PROJECTION AND NEW POINT COORDINATES
        #       - PRODUCES LOCATION OF NEW POINTS
        # CSDL OPERATIONS:
        # create_input: need to use this for setting the CSDL objects of the parameters  (immutable)
        # create_output: explicitly computed output (will use this for intermediate calculations)
        # 
        # register_output: we will use this to register the final set of edge node coordinates
        # new_point_location = matvec(ffd_param, ffd_deltas) + original_points

        self.shape_parameter_model    = ShapeParameterModel(
            shape_parameter_list_input=self.shape_parameter_list,
            shape_parameter_index_input=self.parameter_index_list,
            shape_parametrization=self.shape_param_sps_mat,
        )
        self.mesh_models = []
        for i in range(len(self.rotation_angles)):
            self.mesh_models.append(
                EdgeUpdateModel(
                    ffd_parametrization     = self.ffd_face_sps_mat_list[i],
                    edge_parametrization    = self.edge_param_sps_mat_list[i],
                    initial_edge_coords     = self.initial_edge_coords_instances[i]
                )
            )

        return self.mesh_models, self.shape_parameter_model

    def return_ffd_parameters(self, coordinate_system='cartesian'):
        # OLD METHOD
        # self.initial_edge_coords_instances = []
        # for i in range(len(self.rotation_angles)):
        #     initial_edge_coords_instance = self.get_ffd_edge_old_coords(
        #         output_type=coordinate_system, 
        #         instance=i
        #     )
        #     self.initial_edge_coords_instances.append(initial_edge_coords_instance[:-2]) # TRIMMING OUT ORIGIN

        # NEW METHOD FOR EDGE COORDINATES:
        instances = len(self.rotation_angles)
        self.initial_edge_coords_instances = self.get_original_edge_coords(
            instances=instances
        ) # OUTPUT TYPE IGNORED HERE
        self.ffd_param_dict = {}

        self.ffd_param_dict['shape_parameter_list_input'] = self.shape_parameter_list
        self.ffd_param_dict['shape_parameter_index_input'] = self.parameter_index_list
        self.ffd_param_dict['shape_parametrization'] = self.shape_param_sps_mat
        self.ffd_param_dict['ffd_parametrization'] = self.ffd_face_sps_mat_list
        self.ffd_param_dict['edge_parametrization'] = self.edge_param_sps_mat_list
        self.ffd_param_dict['initial_edge_coordinates'] = self.initial_edge_coords_instances

        '''
        COMMENTS ON FFD PARAMETER DICTIONARY:
        - shape_parameter_list_input: list of shape parameter names (need for csdl.declare_variable)
        - shape_parameter_index_input: list holding info for constant (1) or linear (2) shape parameter
        - shape_parametrization: sparse matrix applied to __ to get the ffd control point deltas
        - THE BELOW QUANTITIES ALL HAVE "n" INSTANCES BASED ON THE LENGTH OF self.rotation_angles INPUT
        - ffd_parametrization: sparse matrix applied to __ to get deltas of mesh vertices
          - these get added to the original mesh vertices
        - edge_parametrization: sparse matrix applied to updated location of mesh vertices to get new edge node locations
        - ffd_cps: initial ffd face control point coordinates
        - initial_edge_coordinates: initial location of all edge nodes
        '''

        return self.ffd_param_dict
   
    # --------------------- MISCELLANEOUS ---------------------
    def extract_point_node_coords(self, points):
         # VECTORIZE &  REMOVE Z-COORDINATE
        dim = 2
        tol = 1e-6
        extracted_point_node_coords     = []

        gmsh_point_node_coords          = self.gmsh_point_nodes[1] # VECTORIZED LIKE [X0 Y0 Z0 X1 Y1 Z1 ...]
        for point in points:
            point_coordinate    = point.return_coordinates()[:dim]
            extracted_point_node_coords.append(point_coordinate)
        
        if self.point_node_coords_to_file:
            file_data = np.array(extracted_point_node_coords).reshape(dim*len(points),1)
            np.savetxt(self.point_node_coords_to_file + '_{}.txt'.format(int(self.instance+1)), file_data)
        
        return extracted_point_node_coords

    def extract_point_node_indices(self, points): # called during GMSH operations in assemble_mesh
        tol = 1e-6
        extracted_point_node_indices    = []

        gmsh_point_node_indices         = self.gmsh_point_nodes[0]
        gmsh_point_node_coords          = self.gmsh_point_nodes[1] # VECTORIZED LIKE [X0 Y0 Z0 X1 Y1 Z1 ...]

        point_coordinates               = self.extract_point_node_coords(points)
        
        for point_coordinate in point_coordinates:
            # print('point_coordinate: ', point_coordinate)
            # implement search method via point coordinates and np.linalg.norm()
            # search within gmsh.model.mesh.getNodes
            for i, node in enumerate(gmsh_point_node_indices):
                diff   = np.array(point_coordinate) - np.array(gmsh_point_node_coords[3*i:3*i+3])
                if np.linalg.norm(diff) < tol:
                    extracted_point_node_indices.append(node)
                    break
        
        if self.point_node_indices_to_file: # if this is True, self.point_node_indices_to_file is a name-string for the file
            np.savetxt(self.point_node_indices_to_file + '.txt', extracted_point_node_indices, fmt='%3i')

        return extracted_point_node_indices

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
                ]
        # self.mesh_nodes  = self.mesh_nodes.reshape((num_pts * 2,))
        return self.mesh_nodes

    def get_node_coordinates(self, points, print_to_file=False):
        if print_to_file:
            self.point_node_coords_to_file = print_to_file

        if not all([isinstance(point, Point) for point in points]):
            raise TypeError('Variables must all be points.')
        for point in points:
            self.points_for_coord_extraction.append(point)

    def get_node_indices(self, points, print_to_file=False): # print_to_file is the file name without extension
        if print_to_file:
            self.point_node_indices_to_file = print_to_file

        if not all([isinstance(point, Point) for point in points]):
            raise TypeError('Variables must all be points.')
        for point in points:
            self.points_for_node_extraction.append(point)

    def reorder_points_to_gmsh(self, coordinate_system='cartesian'):
        # need to reorder points and edges to align with the final result from gmsh
        # POINTS
        # for points use self.point_coordinates and reorder them
        # print(self.point_coordinates)
        gmsh_order_point_coords = []
        gmsh_order_point_node_ind = []
        move_origin, origin_ind = 0, 0
        for point in self.gmsh_point_entities:
            point_info = gmsh.model.mesh.getNodes(
                    dim=point[0], 
                    tag=point[1], 
                    includeBoundary=True, 
                    returnParametricCoord=False
                )
            [point_ind, point_coord] = point_info[:2]

            if np.linalg.norm(point_coord - np.array([0., 0., 0.])) > 1.e-8:
                gmsh_order_point_coords.append(list(point_coord))
                gmsh_order_point_node_ind.append(point_ind[0])
            else:
                move_origin = 1
                origin_ind = point_ind
        
        if move_origin:
            gmsh_order_point_coords.append([0., 0., 0.])
            gmsh_order_point_node_ind.append(origin_ind[0])

        gmsh_order_point_coords = np.array(gmsh_order_point_coords)
        gmsh_order_point_node_ind = np.array(gmsh_order_point_node_ind)

        if coordinate_system == 'polar':
            for i, coords in enumerate(gmsh_order_point_coords):
                gmsh_order_point_coords[i] = [
                    np.arctan2(coords[1], coords[0]),
                    np.linalg.norm(coords[:2]),
                    0.0,
                ]

        return gmsh_order_point_coords, gmsh_order_point_node_ind
        # we can return the points in the order that GMSH now defines them, in either polar or cartesian
        # need this to relate the point objects to the GMSH defined points

    # Ru: seperate the functions for `old_edge_coords` and `edge_deltas` to 
    # make the reusage of the second function more efficient
    def get_ffd_edge_old_coords(self, output_type, instance=0):
        # ordered_coords = self.gmsh_order_point_coords_polar
        ordered_coords = self.mesh_points_instances[instance]
        # old_edge_coords = self.edge_param_sps_mat.dot(ordered_coords)
        old_edge_coords = self.edge_param_sps_mat_list[instance].dot(ordered_coords)
        if output_type == 'cartesian':
            for i in range(int(len(old_edge_coords)/2)):
                old_edge_coords[2*i:2*i+2] = [
                    old_edge_coords[2*i+1] * np.cos(old_edge_coords[2*i]),
                    old_edge_coords[2*i+1] * np.sin(old_edge_coords[2*i]),
                ]
        return old_edge_coords

    def get_original_edge_coords(self, output_type=None, instances=1):
        '''
        DATA NEEDED:
            - self.edge_node_indices
            - self.edge_nodes_instances (length = number of rotation instances)
        '''
        num_edge_nodes = max([max(node) for node in self.edge_node_indices])
        indices = []
        num_edge_groups = len(self.edge_node_indices)
        for edge_group in self.edge_node_indices:
            indices.extend(edge_group)
        indices = np.array(indices)
        
        extraction_indices = []
        for i in range(1, int(num_edge_nodes+1)):
            first_index = np.where(indices == i)[0][0] # if error occurs, check this double indexing
            extraction_indices.append(first_index)
        
        initial_edge_coords_instances = []
        for i in range(instances):
            # FORMATTING COORDINATES IN AN ORDER SIMILAR TO THE INDICES ABOVE
            point_coords = []
            
            for point_coord in self.edge_nodes_instances[i]:
                point_coords.extend(point_coord)
            
            extracted_initial_edge_coords = []
            for j in range(len(extraction_indices)):
                extraction_index = extraction_indices[j]
                extracted_initial_edge_coords.append(point_coords[extraction_index][:2])

            # initial_edge_coords_instance = np.zeros((num_edge_nodes, 2))
            initial_edge_coords_instance = np.reshape(
                np.array(extracted_initial_edge_coords),
                (int(num_edge_nodes*2),)
            )

            # the current approach is in the coordinate system set at the top
            # i.e. initial_edge_coords_instance is in whatever system set above
            if False: 
                if output_type == 'polar': # NEED TO FIX
                    for j in range(num_edge_nodes):
                        coords = initial_edge_coords_instance[j,:]
                        initial_edge_coords_instance[j,:] = [
                            np.arctan2(coords[1], coords[0]),
                            np.linalg.norm(coords)
                        ]

            initial_edge_coords_instances.append(initial_edge_coords_instance)
            # NOTE: may need to revisit to look at the trimming of the origin
        
        return initial_edge_coords_instances

        
    def test_ffd_edge_parametrization_polar(self, delta, output_type, instance=0):
        # if delta.shape != (int(4 * self.num_ffd_faces), 2):
        #     raise TypeError('Shape of delta numpy array must be 4 * # FFD faces * 2')
        delta = delta.reshape((2 * delta.shape[0], ))

        # --- FFD TEST
        # ordered_coords = self.gmsh_order_point_coords_polar
        ordered_coords = self.mesh_points_instances[instance]
        # print(self.ffd_face_sps_mat.shape)
        point_delta = self.ffd_face_sps_mat_list[instance].dot(delta)
        new_coords = ordered_coords + point_delta

        # --- EDGE COORDINATE TEST
        new_edge_coords = self.edge_param_sps_mat_list[instance].dot(new_coords)
        
        if output_type == 'cartesian':
            for i in range(int(len(new_edge_coords)/2)):
                new_edge_coords[2*i:2*i+2] = [
                    new_edge_coords[2*i+1] * np.cos(new_edge_coords[2*i]),
                    new_edge_coords[2*i+1] * np.sin(new_edge_coords[2*i]),
                ]
        old_edge_coords = self.get_ffd_edge_old_coords(output_type=output_type, instance=instance)
        edge_deltas = new_edge_coords - old_edge_coords
        edge_indices = []
        for i in range(len(edge_deltas)):
            edge_indices.append(i)

        return edge_deltas # need to convert to cartesian

