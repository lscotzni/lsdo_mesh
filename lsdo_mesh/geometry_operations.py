
import numpy as np
from lsdo_mesh.mesh_ffd_objects import *

# TO-DO: instead of input_type, check the object type in entities
def rotate_old(entities=None, angle=[], center_point = None, copy = True, input_type = 'points'):
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
            p0  = np.array(list(vars(point)['properties'].values())[:3])

            rot_matrix = np.array([np.cos(angle[2]), -np.sin(angle[2]), np.sin(angle[2]), np.cos(angle[2])]).reshape((2,2))

            p1  = [np.dot(rot_matrix[0,:], p0[:2]), np.dot(rot_matrix[1,:], p0[:2]), cp[2]]
            # print(p0)
            # print(p1)
            # exit()

            if copy == True:
                rotated_points.append(Point(p1[0], p1[1], ms=ms))
            elif copy == False:
                pass
            
        return rotated_points

def rotate(entities=None, angle=[], center_point = None, copy = True):
    
    if entities is None:
        raise KeyError('No objects to rotate.')

    first_entity = entities[0]
    first_entity_type = type(entities[0])
    # can also check using if not isinstance(entities[0], entity)
    for entity in entities[1:]:
        if not isinstance(first_entity, type(entity)):
            raise KeyError('Inconsistent geometry types for rotation.')
    
    if center_point == None: 
        cp      = [0., 0., 0.] # DEFAULT TO ORIGIN
    else:
        cp      = list(vars(center_point)['properties'].values())[:3]

    if isinstance(first_entity, Point) or isinstance(first_entity, Vertex):
    # if first_entity_type is Point or first_entity_type is Vertex:
        # geometry_type = 
        rotated_points = []
        for point in entities:
            if isinstance(point, Point):
                ms  = list(vars(point)['properties'].values())[3]
            p0  = np.array(list(vars(point)['properties'].values())[:3])

            rot_matrix = np.array([np.cos(angle[2]), -np.sin(angle[2]), np.sin(angle[2]), np.cos(angle[2])]).reshape((2,2))

            # p1  = [np.dot(rot_matrix[0,:], p0[:2]) - cp[0], np.dot(rot_matrix[1,:], p0[:2])  - cp[1], 0.] + cp
            p1 = [
                np.dot(rot_matrix[0,:], p0[:2]) - cp[0],
                np.dot(rot_matrix[1,:], p0[:2]) - cp[1],
                0.
            ] + cp

            if copy == True:
                if isinstance(first_entity, Point):
                    rotated_points.append(Point(p1[0], p1[1], ms=ms))
                elif isinstance(first_entity, Vertex):
                    rotated_points.append(Vertex(p1[0], p1[1]))
            elif copy == False:
                pass
            
        return rotated_points
    # elif isinstance(first_entity, Surface) or isinstance(first_entity, Face):
