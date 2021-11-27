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

    def initialize(self, *args, ms = 0.1, physical_group = False):
        props = self.properties
        props['x1'] = args[0]
        props['x2'] = args[1]
        props['x3'] = 0.
        props['ms'] = ms

    def assemble_self(self):
        return self.mesh.add_point(self.properties['x1'], self.properties['x2'], self.properties['x3'], self.properties['ms'])

class Curve(Entity):

    def initialize(self, *args, curve_type='line', coord_sys='rectilinear', physical_group=False):
        if curve_type != 'line' and curve_type != 'arc':
            raise KeyError('Line types besides arc and line have not been set up.')
        # error checking needed for later
        self.children = list(args) 
        props = self.properties

        if coord_sys == 'rectilinear':
            props['coord_sys'] = 0
            if curve_type == 'line':
                props['type'] = 0
            elif curve_type == 'arc':
                props['type'] = 1
        elif coord_sys == 'polar':
            props['coord_sys'] = 1
            start_point = list(vars(self.children[0])['properties'].values())
            end_point = list(vars(self.children[1])['properties'].values())
            # [r, theta, z, ms]
            
            if start_point[0] == end_point[0] and start_point[1] != end_point[1]:
                # curve of constant radius (arc, requires 3 points)
                if len(self.children) != 3:
                    raise KeyError('Missing center point for circle arc.')
                props['type'] = 1

            elif start_point[0] != end_point[0] and start_point[1] == end_point[1]:
                # curve of constant theta (line, requires 2 points)
                props['type'] = 0
            elif start_point[0] != end_point[0] and start_point[1] != end_point[1]:
                raise KeyError('Cannot create curve in polar with varying R and THETA.')

    def assemble_self(self):
        return self.mesh.add_curve(
            self.children, 
            curve_type=self.properties['type'], 
            coord_sys=self.properties['coord_sys']
        )

# class BooleanCurve(...)

class Surface(Entity):

    def initialize(self, *args, input_type='curves', physical_group = 0):
        props = self.properties
        if input_type == 'curves':
            props['type'] = 0
            args    = args[0]
        elif input_type == 'polygon':
            curves = []
            for ind in range(len(args[0])):
                if ind < len(args[0]) - 1:
                    curve = Curve(args[0][ind], args[0][ind + 1])
                else:
                    curve = Curve(args[0][ind], args[0][0])
                curves.append(curve)
            args = curves
        else:
            raise KeyError('Only curves and polygon have been implemented.')

        self.children       = list(args)
        self.physical_group = physical_group

    def assemble_self(self):
        return self.mesh.add_surface(self.children, self.physical_group)

    def set_physical_group(self, physical_group=(0,None)):
        pass

class BooleanSurface(Entity):
    
    def initialize(self, objects, tools, operation, removeObject=True, removeTool=True, geometry='surfaces', physical_group = 0):
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

    def initialize(self, *args, input_type='edges'):
        props = self.properties

        if isinstance(args[0], Point):
            p = list(vars(args[0])['properties'].values())

            if input_type == 'cartesian':
                p1 = Vertex(p[0] * (1 + 0.01), p[1] * (1 + 0.01))
                p2 = Vertex(p[0] * (1 - 0.01), p[1] * (1 + 0.01))
                p3 = Vertex(p[0] * (1 - 0.01), p[1] * (1 - 0.01))
                p4 = Vertex(p[0] * (1 + 0.01), p[1] * (1 - 0.01))

                # Face([p1, p2, p3, p4], input_type='polygon')

            elif input_type == 'polar':
                r = np.linalg.norm(p[:2])
                theta = np.arctan2(p[1], p[0])
                p1 = Vertex(r * (1 + 0.01), theta * (1 - 0.01), mode='polar')
                p2 = Vertex(r * (1 + 0.01), theta * (1 + 0.01), mode='polar')
                p3 = Vertex(r * (1 - 0.01), theta * (1 + 0.01), mode='polar')
                p4 = Vertex(r * (1 - 0.01), theta * (1 - 0.01), mode='polar')

                # Face([p1, p2, p3, p4], input_type='polygon')
            vertices = [p1, p2, p3, p4]
            edges = []
            for ind in range(len(vertices)):
                if ind < len(vertices) - 1:
                    edge = Edge(vertices[ind], vertices[ind + 1])
                else:
                    edge = Edge(vertices[ind], vertices[0])
                edges.append(edge)
            args = edges
            
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
        print(self.children)

    def assemble_self(self):
        return self.mesh.define_face(self.children)

    def add_shape_parameter(self, name=None, bounds=None, axis=None, degree=None):
        pass

class Block(Entity):
    def initialize(self,*args):
        pass