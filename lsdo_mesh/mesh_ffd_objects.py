import numpy as np

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

    def initialize(self, *args, ms = 0.1, mode='cartesian', physical_group = False):
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
        self.mode = mode

    def assemble_self(self):
        return self.mesh.add_point(self.properties['x'], self.properties['y'], self.properties['z'], self.properties['ms'])
    
    def return_coordinates(self, output_type='cartesian'):
        if output_type is 'cartesian':
            return list(self.properties.values())[:3]
        elif output_type is 'polar':
            cartesian_coords = list(self.properties.values())[:3]
            polar_coords = [
                np.arctan2(cartesian_coords[1], cartesian_coords[0]),
                np.linalg.norm(cartesian_coords[:2]),
                cartesian_coords[2]
            ]
            return polar_coords[:3]

class Curve(Entity):

    def initialize(self, *args, curve_type='line', physical_group=False):
        if curve_type != 'line' and curve_type != 'arc':
            raise KeyError('Line types besides arc and line have not been set up.')
        # error checking needed for later
        self.children = list(args) 

        props = self.properties
        if curve_type == 'line':
            props['type'] = 0
        elif curve_type == 'arc':
            props['type'] = 1

        self.physical_group = physical_group

    def assemble_self(self):
        return self.mesh.add_curve(self.children, curve_type=self.properties['type'], physical_group=self.physical_group)

# class BooleanCurve(...)

class Surface(Entity):

    def initialize(self, *args, input_type='curves', physical_group=False):
        props = self.properties
        curve_loop_lengths = []
        if input_type == 'curves':
            props['type'] = 0
            args    = args[0]
            curve_loop_lengths.append(len(args))
            

        elif input_type == 'polygon':
            curves = []
            for ind in range(len(args[0])):
                if ind < len(args[0]) - 1:
                    curve = Curve(args[0][ind], args[0][ind + 1])
                else:
                    curve = Curve(args[0][ind], args[0][0])
                curves.append(curve)
            args = curves
            curve_loop_lengths.append(len(args))
        elif input_type == 'curve_loops': # IMPLIES AT LEAST 2 DISTINCT SETS OF CURVES
            props['type'] = 1
            
            curves = []
            for loop in args: # args formatted as ([[...]], [[...]], ...)
                curve_loop_lengths.append(len(loop[0]))
                for curve in loop[0]: # loop formatted as [[...]]; loop[0] is the internal list
                    curves.append(curve)
            args = curves
            # turns embedded lists in original args to single list of curves
            # need a counter to signify which curves exist for which curve loop

        else:
            raise KeyError('Only curves, polygon and curve loops have been implemented.')

        self.children       = list(args)
        self.physical_group = physical_group
        # the format for a physical group is (ind, name) or ind; 
        self.curve_loops = curve_loop_lengths

    def assemble_self(self):
        return self.mesh.add_surface(self.children, self.curve_loop_lengths, self.physical_group)

class BooleanSurface(Entity):
    
    def initialize(self, objects, tools, operation, removeObject=True, removeTool=True, geometry='surfaces', physical_group=False):
        props           = self.properties
        self.children   = []
        self.children.extend(objects)
        self.children.extend(tools)

        props['objects']        = len(objects)
        props['tools']          = len(tools)
        props['operation']      = operation
        props['removeObject']   = removeObject
        props['removeTool']     = removeTool
        if isinstance(objects[0], Surface):
            props['dimension'] = 2
        else:
            raise KeyError('Boolean operations currently defined for only surfaces.')

        self.physical_group     = physical_group

    def assemble_self(self):
        return self.mesh.add_boolean(self.children, tuple(self.properties.values()), self.physical_group)

class Volume(Entity):
    def initialize(self,*args):
        pass

# FFD

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
    
    def return_coordinates(self):
        return list(self.properties.values())

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
        return self.mesh.define_edge(self.children, edge_type = self.properties['type'])

class Face(Entity):

    def initialize(self, *args, input_type='edges', delta=1.):
        props = self.properties
        self.embedded_points = []

        if isinstance(args[0], Point) or all([isinstance(arg, Point) for arg in args[0]]):
            delta /= 100.
            if isinstance(args[0], Point):
                self.embedded_points.append(args[0])
                delta /= 100.
                p = list(vars(args[0])['properties'].values())

                if input_type == 'cartesian':
                    p1 = Vertex(p[0] * (1 + delta), p[1] * (1 + delta))
                    p2 = Vertex(p[0] * (1 - delta), p[1] * (1 + delta))
                    p3 = Vertex(p[0] * (1 - delta), p[1] * (1 - delta))
                    p4 = Vertex(p[0] * (1 + delta), p[1] * (1 - delta))
                    # print('face worked')
                    vertices = [p1, p2, p3, p4]
                    self.face_vertices = vertices
                    Face(vertices, input_type='quad')
                    
                elif input_type == 'polar':
                    r = np.linalg.norm(p[:2])
                    theta = np.arctan2(p[1], p[0])
                    p1 = Vertex(r * (1 + delta), theta * (1 - delta) - 0.1, mode='polar')
                    p2 = Vertex(r * (1 + delta), theta * (1 + delta) + 0.1, mode='polar')
                    p3 = Vertex(r * (1 - delta), theta * (1 + delta) + 0.1, mode='polar')
                    p4 = Vertex(r * (1 - delta), theta * (1 - delta) - 0.1, mode='polar')

                    vertices = [p1, p2, p3, p4]
                    self.face_vertices = vertices
                # ADD SOMETHING THAT SAYS ARGS = VERTICES
                # COULD BE AS EASY AS args = vertices
                # the Point method to create Faces works well but if we 
                # define an ARBITRARY face we need to be able to extract
                # the embedded points and add them to a separate variable
                # like self.embedded_points
            # elif len(args[0]) != 1:
            else:
                self.embedded_points.extend(args[0])
                coords = []
                np.array(coords.extend([arg.return_coordinates()[:2] for arg in args[0]]))
                # line above will return CARTESIAN COORDINATES for simplicity (do not change)
                if input_type == 'cartesian':
                    # these two lines below are wrong
                    min_pt = [min(coords[:,0]), min(coords[:,1])]
                    max_pt = [max(coords[:,0]), max(coords[:,1])]
                    
                    p1 = Vertex(max_pt[0] * (1 + delta), max_pt[1] * (1 + delta))
                    p2 = Vertex(min_pt[0] * (1 - delta), max_pt[1] * (1 + delta))
                    p3 = Vertex(min_pt[0] * (1 - delta), min_pt[1] * (1 - delta))
                    p4 = Vertex(max_pt[0] * (1 + delta), min_pt[1] * (1 - delta))
                    vertices = [p1, p2, p3, p4]
                    self.face_vertices = vertices
                    Face(vertices, input_type='quad')

                elif input_type == 'polar':
                    max_pt = [max([np.linalg.norm(coord) for coord in coords]), max([np.arctan2(coord[1], coord[0]) for coord in coords])]
                    min_pt = [min([np.linalg.norm(coord) for coord in coords]), min([np.arctan2(coord[1], coord[0]) for coord in coords])]
                
                    p1 = Vertex(max_pt[0] * (1 + delta), max_pt[1] * (1 + delta), mode='polar')
                    p2 = Vertex(max_pt[0] * (1 + delta), min_pt[1] * (1 - delta), mode='polar')
                    p3 = Vertex(min_pt[0] * (1 - delta), min_pt[1] * (1 - delta), mode='polar')
                    p4 = Vertex(min_pt[0] * (1 - delta), max_pt[1] * (1 + delta), mode='polar')

                    vertices = [p1, p2, p3, p4]
                    self.face_vertices = vertices
            args = vertices

        elif input_type == 'edges':
            props['type'] = 0
            args = args[0]
        elif input_type == 'quad': # general quad
            edges = []
            for ind in range(len(args[0])):
                if ind < len(args[0]) - 1:
                    edge = Edge(args[0][ind], args[0][ind + 1])
                else:
                    edge = Edge(args[0][ind], args[0][0])
                edges.append(edge)
            args = edges
        elif input_type == 'circle_sector':
            pass
        else:
            raise KeyError('Only edges and polygon have been implemented.')

        self.children = list(args)

    def assemble_self(self):
        return self.mesh.define_face(self.children)

    def add_shape_parameter(self, name=None, axis=None, def_type=None):
        self.parameters = {
            'name': name,
            'axis': axis,
            'type': def_type
        }
        # need to get parametrization of point based on external vertices
        pts = np.array(self.return_coordinates())
        # we can take the 4 vertices of the face and find the parametric coordinates of
        # the point of interest so that any movement of those 4 vertices will be seen in
        # the node of interest
        if axis is 'x' or axis is 'y':
            pass # working on polar first
        elif axis is 'r' or axis is 'theta': # theta for u, r for v
            P00 = np.min(pts, axis=0)
            P11 = np.max(pts, axis=0)
            P01 = [np.min(pts[:,0]), np.max(pts[:,1])]
            P10 = [np.max(pts[:,0]), np.min(pts[:,1])]
            # print(P00, P11, P01, P10)

        return self.parameters

    def add_embedded_points(self, points):
        # this method is so that we can arbitrarily add points
        # as properties of this FFD face
        # add check to make sure that the coordinates lie within the face
        for point in points:
            # if point in range of vertices:
                # this if clause can check if x & y of point are within
                # the range of the vertex x & y points 
                # self.embedded_points.append(point)
            pass

        pass
        
    def return_coordinates(self, coordinate_system='cartesian'):
        vertex_coords = []
        for vertex in self.face_vertices:
            vertex_coords.append(vertex.return_coordinates()[:2])
            # print(vertex.return_coordinates()[:2])
            # THIS RETURNS COORDINATES IN X-Y SPACE
        if coordinate_system is 'polar':
            for i, vertex in enumerate(vertex_coords):
                vertex_coords[i] = [
                    np.arctan2(vertex[1], vertex[0]),
                    np.linalg.norm(vertex)
                ]
        return vertex_coords

class Block(Entity):
    def initialize(self,*args):
        pass