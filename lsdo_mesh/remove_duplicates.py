import numpy as np
import time

'''
-------------------- REMOVING DUPLICATE POINTS --------------------
'''

def remove_duplicate_points(points=None, point_mesh_size=None, point_rotate_instance=None, curves=None, eps=1.e-8):
    dummy_mesh_size     = False

    # -- REMOVE --
    # This block is how to deal with code for FFD duplicates (since we don't have a mesh size)
    if point_mesh_size is None:
        point_mesh_size = np.zeros(len(points))
        dummy_mesh_size = True

    t1 = time.time()
    point_indices_to_connect = []
    iterated_indices = []
    for i in range(len(points)):
        if i in iterated_indices:
            continue
        for j in range(i+1,len(points)):
            if np.linalg.norm(np.array(points[i,:] - points[j,:])) < eps:
                point_indices_to_connect.append((i,j))
                iterated_indices.append(j)
                # if j not in point_indices_to_connect:
                #     point_indices_to_connect.append((i,j))

    t2 = time.time()
    print('duplicate indices length', len(point_indices_to_connect))
    print('point connection time:', t2 - t1)

    # NOTE: full motor takes 292 seconds and point_indices_to_connect is 178760 long

    # print(point_indices_to_connect[:10])
    print('Completed search of point indices connections.')
    # print('{} points were connected out of {} points.'.format(len(point_indices_to_connect), len(points)))
    
    new_points, new_point_mesh_size, new_point_rotate_instances = combine_points(points, point_mesh_size, point_rotate_instance, point_indices_to_connect)
    t3 = time.time()
    print('Points combined.')
    print('num new points', new_points.shape[0])
    print('duplicate removal time:', t3 - t2)
    # exit()

    new_curves = np.zeros(len(curves), dtype=int)
    for i,new_point in enumerate(new_points):
        for j,point in enumerate(points):
            if np.linalg.norm(point - new_point) < eps:
                new_curves[j] = curves[i]
    print('New curves from point removal created.')

    if dummy_mesh_size == False:
        return new_points, new_point_mesh_size, new_point_rotate_instances, new_curves
    elif dummy_mesh_size == True:
        return new_points, new_curves

def combine_points(points, point_mesh_size, point_rotate_instance, point_indices_to_connect):
    new_points = np.delete(points,[indices[1] for indices in point_indices_to_connect],axis=0)
    new_point_mesh_size = np.delete(point_mesh_size,[indices[1] for indices in point_indices_to_connect],axis=0)
    new_point_rotate_instances = np.delete(point_rotate_instance,[indices[1] for indices in point_indices_to_connect],axis=0)
    return new_points, new_point_mesh_size, new_point_rotate_instances



'''
-------------------- REMOVING DUPLICATE CURVES --------------------
'''
def remove_duplicate_curves(curves=None, curve_indices=None, curve_type=None, curve_physical_groups=None, surfaces=None):
    dummy_curve_type     = False

    # This block is how to deal with code for FFD duplicates (since we don't have a curve type)
    if curve_type is None:
        curve_type = np.zeros(len(surfaces))
        dummy_curve_type = True

    curve_indices_to_connect = []
    for i in range(len(curve_indices)):
        for j in range(i+1,len(curve_indices)):
            if curves[curve_indices[i,0]]==curves[curve_indices[j,0]] and curves[curve_indices[i,0]+1]==curves[curve_indices[j,0]+1] \
            and curve_type[i] == curve_type[j]:
                curve_indices_to_connect.append((i,j))
            if curves[curve_indices[i,0]]==curves[curve_indices[j,0]+1] and curves[curve_indices[i,0]+1]==curves[curve_indices[j,0]] \
             and curve_type[i] == curve_type[j]:
                curve_indices_to_connect.append((i,j))
    # NOTE:
    # THIS curve_indices_to_connect METHOD IS BAD WHEN WE ARE LOOKING FOR THE LAST INDEX; THIS DIFFERS
    # FOR CIRCLE ARCS BC THE LAST POINT IS THE CENTER POINT, BUT FOR ALL OTHER CURVES LIKE LINES AND 
    # SPLINES THIS IS NOT THE CASE -----> FIX

    new_curves, new_curve_indices, new_curve_types, new_curve_physical_groups = combine_curves(
        curves,curve_indices, curve_type, curve_indices_to_connect, curve_physical_groups
        )

    new_surfaces    = np.zeros((len(surfaces),), dtype = int) 
    new_surfaces[:] = surfaces

    for i, index in enumerate(curve_indices_to_connect):
        remove_curve                    = np.where(surfaces == index[1])
        new_surfaces[remove_curve]      = index[0]

        # if curve_physical_groups[index[0]]:
        #     curve_physical_groups[index[1]] = False
        # elif curve_physical_groups[index[1]]:
        #     curve_physical_groups[index[0]] = False

    unique_curve_indices_connect        = sorted(list(set([index[0] for index in curve_indices_to_connect])))
    temp_unique_curve_indices_connect   = np.array([i for i in unique_curve_indices_connect])
    unique_surfaces                     = list(set(new_surfaces))

    if dummy_curve_type == False:
        return new_curves, new_curve_indices, new_curve_types, new_curve_physical_groups, new_surfaces, unique_surfaces
    elif dummy_curve_type == True:
        return new_curves, new_curve_indices, new_surfaces, unique_surfaces

def combine_curves(curves, curve_indices,curve_type, curve_indices_to_connect, curve_physical_groups):
    # NEW CURVES
    index_list_unformatted = [list(np.arange(curve_indices[indices[1]][0],curve_indices[indices[1]][1],dtype = int)) for indices in curve_indices_to_connect]
    index_list = [ind for inds_sublist in index_list_unformatted for ind in inds_sublist] # converting from lists within list to a single list of numbers
    new_curves  = np.delete(curves,index_list,axis=0)
    # NEW CURVE INDICES & CURVE TYPES
    temp_new_curve_indices  = np.delete(curve_indices,[indices[1] for indices in curve_indices_to_connect],axis=0) # correct size, wrong indices
    new_curve_types         = np.delete(curve_type,[indices[1] for indices in curve_indices_to_connect],axis=0)
    new_curve_physical_groups = list(np.delete(curve_physical_groups,[indices[1] for indices in curve_indices_to_connect]))
    
    new_curve_indices       = []
    counter                 = 0

    for i, curve_type_int in enumerate(new_curve_types):
        # print(i)
        temp_len            = temp_new_curve_indices[i,1] - temp_new_curve_indices[i,0] # number of points for curves in iterations
        new_curve_indices.extend([counter,counter + temp_len])
        counter             += temp_len

    new_curve_indices       = np.array(new_curve_indices,dtype = int).reshape(len(new_curve_types),2)

    return new_curves, new_curve_indices, new_curve_types, new_curve_physical_groups

'''
NEED TO IMPLEMENT A REMOVE_DUPLICATE_SURFACES
'''