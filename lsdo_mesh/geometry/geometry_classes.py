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

'''
=============================================================================================================================================
------------------------------------------------------------ FFD CLASSES/METHODS ------------------------------------------------------------
=============================================================================================================================================
'''

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