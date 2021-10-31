import gmsh
import numpy as np
import sys

from lsdo_mesh.mesh.remove_duplicates import remove_duplicate_points, remove_duplicate_curves
# from lsdo_mesh.geometry.geometry_classes import Entity, Point, Curve, Surface, Vertex, Edge, 

class Entity(object):

    def __init__(self, *args, **kwargs):
        self.children = []
        self.properties = dict()
        self.mesh = None
        self.initialize(*args, **kwargs)
        self.id = None

    def assemble(self, mesh):
        self.mesh = mesh

        for child in self.children:
            child.assemble(mesh)
            # THIS IS A RECURSIVE LOOP THAT FEEDS INTO THE OBJECTS FOR THE DIMENSIONS BELOW
        
        self.id = self.assemble_self()

class Point(Entity):

    def initialize(self, *args, ms = 0.1, mode='cartesian'):
        props = self.properties
        if mode is 'cartesian':
            props['x'] = args[0]
            props['y'] = args[1]
            props['z'] = 0.
        elif mode is 'polar':
            r = args[0]
            theta = args[1]
            props['x'] = r * np.cos(theta)
            props['y'] = r * np.sin(theta)
            props['z'] = 0.
        props['ms'] = ms

    def assemble_self(self):
        return self.mesh.add_point(self.properties['x'], self.properties['y'], self.properties['z'], self.properties['ms'])

class Curve(Entity):

    def initialize(self, *args, curve_type='line'):
        if curve_type != 'line' and curve_type != 'arc':
            raise KeyError('Line types besides arc and line have not been set up.')
        # error checking needed for later
        self.children = list(args) 

        props = self.properties
        if curve_type == 'line':
            props['type'] = 0
        elif curve_type == 'arc':
            props['type'] = 1

    def assemble_self(self):
        return self.mesh.add_curve(self.children,curve_type = self.properties['type'])

class Surface(Entity):

    def initialize(self, *args, input_type='curves'):
        props = self.properties
        if input_type == 'curves':
            props['type'] = 0
        elif input_type == 'polygon':
            curves = []
            for ind in range(len(args)):
                if ind < len(args) - 1:
                    curve = Curve(args[ind], args[ind + 1])
                else:
                    curve = Curve(args[ind], args[0])
                curves.append(curve)
            args = curves
        else:
            raise KeyError('Only curves and polygon have been implemented.')

        self.children = list(args)

    def assemble_self(self):
        return self.mesh.add_surface(self.children)

class Volume(Entity):
    def initialize(self,*args):
        pass

# class BooleanCurve(...)

class BooleanSurface(Entity):
    
    def initialize(self, objects, tools, operation, removeObject=True, removeTool=True, geometry='surfaces'):
        props           = self.properties
        self.children   = []
        self.children.extend(objects)
        self.children.extend(tools)

        props['objects']        = len(objects)
        props['tools']          = len(tools)
        props['operation']      = operation
        props['removeObject']   = removeObject
        props['removeTool']     = removeTool
        if geometry == 'surfaces':
            props['dimension'] = 2
        else:
            raise KeyError('Boolean operations currently defined for only surfaces.')

    def assemble_self(self):
        return self.mesh.add_boolean(self.children,list(self.properties.values()))

class Mesh(object):

    def __init__(self,name = 'Mesh',popup = False):
        self.name       = name
        self.popup      = popup
        
        self.top_entities           = []

        self.point_coordinates      = []
        self.point_mesh_size        = []

        self.curves                 = []
        self.curve_indices          = []
        self.curve_type             = []

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


    def add_entity(self, entity):
        self.top_entities.append(entity)

    def add_point(self, x, y, z, mesh_size):
        self.point_coordinates.append([x, y, z])
        self.point_mesh_size.append(mesh_size)
        return len(self.point_coordinates) - 1

    def add_curve(self, children, curve_type):

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
        return len(self.curve_type) - 1

    def add_surface(self, children):

        self.surfaces.extend([curve.id for curve in children])

        if self.surface_indices == []:
            self.surface_indices.append([0, len(children)])
        else:
            self.surface_indices.append([
                self.surface_indices[-1][-1],
                self.surface_indices[-1][-1] + len(children),
            ])

        return len(self.surface_indices) - 1

    def add_boolean(self,children,properties):

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
        return len(self.boolean_operations) - 1


    def assemble(self):
        # (1) recursively assemble entities data structures
        # (2) remove duplicates
        # (3) gmsh
        # initialize gmsh
        # # create points, curves, surfaces in order (ignoring surfaces defined as a result of boolean operations)
        # # run the boolean operations in order
        # # finalize gmsh 

        for entity in self.top_entities:
            entity.assemble(self)

        # print(self.boolean_operations)
        # print(self.boolean_entities)
        # print(self.boolean_object_indices)
        # print(self.boolean_tool_indices)
        # print(self.boolean_remove_object)
        # print(self.boolean_tool_object)
        # print(self.boolean_parameters)
        # exit()

        # converting to numpy arrays
        self.point_coordinates      = np.array(self.point_coordinates)
        self.point_mesh_size        = np.array(self.point_mesh_size)

        self.curves                 = np.array(self.curves)
        self.curve_indices          = np.array(self.curve_indices)
        self.curve_type             = np.array(self.curve_type)

        self.surfaces               = np.array(self.surfaces)
        self.surface_indices        = np.array(self.surface_indices)
        self.surface_type           = np.array(self.surface_type)

        self.boolean_operations     = np.array(self.boolean_operations)
        self.boolean_entities       = np.array(self.boolean_entities)
        self.boolean_object_indices = np.array(self.boolean_object_indices)
        self.boolean_tool_indices   = np.array(self.boolean_tool_indices)
        self.boolean_remove_object  = np.array(self.boolean_remove_object)
        self.boolean_tool_object    = np.array(self.boolean_tool_object)
        self.boolean_parameters     = np.array(self.boolean_parameters)


        # add something like above that steps through boolean operations
        # array for type of boolean operation, indices of surfaces, etc

        # removing duplicate points
        new_points, new_point_mesh_size, new_curves_temp = remove_duplicate_points(self.point_coordinates,self.point_mesh_size,self.curves)

        # removing duplicate curves
        new_curves, new_curve_indices, new_curve_types, new_surfaces, unique_surfaces = remove_duplicate_curves(new_curves_temp,self.curve_indices,self.curve_type,self.surfaces)

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal",1)
        gmsh.model.add(self.name) 
        occ_kernel      = gmsh.model.occ

        # CREATE POINTS
        for i, point in enumerate(new_points):
            occ_kernel.addPoint(point[0], point[1], point[2], new_point_mesh_size[i], i + 1)

        # boolean operations for points

        # CREATE CURVES
        print(unique_surfaces)
        for i, curve in enumerate(new_curve_indices):
            if new_curve_types[i] == 0: # LINE
                print('i', i, 'Line')
                occ_kernel.addLine(new_curves[curve[0]] + 1, new_curves[curve[1] - 1] + 1, unique_surfaces[i] + 1)

            elif new_curve_types[i] == 1: # CIRCLE ARC
                print('i', i, 'Arc')
                occ_kernel.addCircleArc(new_curves[curve[0]] + 1, new_curves[curve[1] - 1] + 1, new_curves[curve[1] - 2] + 1, unique_surfaces[i] + 1)
            else:
                pass

        #  boolean operations for curves
        
        # CREATE SURFACES
        for i, surface in enumerate(self.surface_indices):
            
            curve_input = list(new_surfaces[np.arange(surface[0],surface[1])]+1)
            curveloop = occ_kernel.addCurveLoop(curve_input)  # fix the  ccc via Geometry.OCCAutoFix = 0 later
            occ_kernel.addPlaneSurface([curveloop],i+1)

        # EXECUTE BOOLEAN OPERATIONS
        for i, parameters in enumerate(self.boolean_parameters):
            print(parameters)
            print(i)
            if parameters[0] == 'subtract':
                # occ_kernel.cut(,removeObject=parameters[2],removeTool=parameters[3]) 
                occ_kernel.cut(
                    [(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_object_indices[i][0],self.boolean_object_indices[i][1])],
                    [(parameters[3],self.boolean_entities[j]+1) for j in np.arange(self.boolean_tool_indices[i][0],self.boolean_tool_indices[i][1])],
                    removeObject=parameters[1],
                    removeTool=parameters[2]
                )
                pass
            else:
                raise KeyError('operation has not been integrated yet')

            # motor.cut([(2,1)],rotor_slot_bool,5 + 3*p + 2*s + 2,removeObject=True,removeTool=False)



        # exit()
        occ_kernel.synchronize()

        gmsh.model.mesh.generate(2)
        gmsh.write(self.name +  '.msh')
        # if '-nopopup' not in sys.argv:
        if self.popup == True:
            gmsh.fltk.run()
        gmsh.finalize()

    def boolean(self,objects, tools, operation, removeObject=True, removeTool=True):
        # fuse, intersect, cut and fragment take all of the same options:
        # object tag, tool tag, remove object (T/F), remove tool (T/F)
        
        self.boolean_type.append(operation)
        self.boolean_objects.extend([entity.id for entity in objects])
        # self.boolean_object_indices()

        self.boolean_tools.extend([entity.id for entity in tools])
        # self.boolean_tool_indices()
        # WE NEED TO FEED INTO add_entity FOR IT TO BE COMBINED
        return len(self.boolean_type) - 1



    def subtract(self,mainObject = [], toolObjects = []):
        print(vars(mainObject[0]))
        print(vars(toolObjects[0]))
        # exit()
        if self.boolean_surfaces == []:
            self.boolean_surfaces.append()
        else:
            pass

'''
-------------------- GEOMETRICAL OPERATIONS --------------------
'''

def rotate(entities=None, angle=[], center_point = None, copy = True, input_type = 'points'):
    if entities is None:
        entities =  []
    if center_point == None: 
        cp      = [0., 0., 0.] # DEFAULT TO ORIGIN
    else:
        cp      = list(vars(center_point)['properties'].values())[:3]
    
    if input_type == 'points':

        rotated_points = []
        for point in entities:
            ms  = list(vars(point)['properties'].values())[3]
            p0  = list(vars(point)['properties'].values())[:3]
            r   = np.sqrt(np.sum([(p0[i] - cp[i])**2 for i in range(len(cp))]))
            base_angle  = np.arctan2(p0[1],p0[0])

            p1  = [cp[0] + r * np.cos(angle[2] + base_angle), cp[1] + r * np.sin(angle[2] + base_angle), cp[2]]

            if copy == True:
                rotated_points.append(Point(p1[0], p1[1], ms=ms))
            elif copy == False:
                pass
            
        return rotated_points

'''
-------------------- BOOLEAN OPERATIONS --------------------
'''
def subtract(mainObject = [], toolObjects = []):
    print(vars(mainObject[0]))
    print(vars(toolObjects[0]))
    # exit()
    # if self.boolean_surfaces == []:
    #     self.boolean_surfaces.append()
    # else:
    #     pass

'''
=============================================================================================================================================
------------------------------------------------------------ FFD CLASSES/METHODS ------------------------------------------------------------
=============================================================================================================================================
'''

# class FFD(object):
#     ...
# inputting design variables to determine shift or stretch
#

class Vertex(Entity):

    def initialize(self, *args, mode='cartesian'):
        props = self.properties
        if mode is 'cartesian':
            props['x'] = args[0]
            props['y'] = args[1]
            props['z'] = 0.
        elif mode is 'polar':
            r = args[0]
            theta = args[1]
            props['x'] = r * np.cos(theta)
            props['y'] = r * np.sin(theta)
            props['z'] = 0.

    def assemble_self(self):
        return self.mesh.define_vertex(self.properties['x'], self.properties['y'], self.properties['z'])

class Edge(Entity):

    def initialize(self, *args, edge_type='line'):
        if edge_type != 'line' and edge_type != 'arc':
            raise KeyError('Line types besides arc and line have not been set up.')
        # error checking needed for later
        self.children = list(args) 

        props = self.properties
        if edge_type == 'line':
            props['type'] = 0
        elif edge_type == 'arc':
            props['type'] = 1

    def assemble_self(self):
        return self.mesh.define_edge(self.children,edge_type = self.properties['type'])

class Face(Entity):

    def initialize(self, *args, input_type='edges'):
        props = self.properties
        if input_type == 'edges':
            props['type'] = 0
        elif input_type == 'polygon':
            edges = []
            for ind in range(len(args)):
                if ind < len(args) - 1:
                    edge = Edge(args[ind], args[ind + 1])
                else:
                    edge = Edge(args[ind], args[0])
                edges.append(edge)
            args = edges
        else:
            raise KeyError('Only edges and polygon have been implemented.')

        self.children = list(args)

    def assemble_self(self):
        return self.mesh.define_face(self.children)

class Block(Entity):
    def initialize(self,*args):
        pass

class FFD(Entity):
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

    def add_design_variable(self):
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

        print('---')
        print('vertices:', self.vertices)
        print('edges:', self.edges)
        print('edge indices:', self.edge_indices)
        print('faces:', self.faces)
        print('face indices:', self.face_indices)
        print('---')

        ffd_grid    = self.create_ffd_grid(self.ffd_grid_steps)


