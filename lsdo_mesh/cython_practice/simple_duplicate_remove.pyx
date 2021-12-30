import numpy
cimport numpy
import time

cdef numpy.ndarray full_list, point_indices_to_connect

cdef float eps, t1, t2, t3, t_connect, t_remove

eps  = 1.e-8

def simple_duplicate_remove(numpy.ndarray full_list):
    cdef list point_indices_to_connect
    cdef list iterated_points
    cdef numpy.ndarray reduced_list
    cdef int full_list_length

    full_list_length = full_list.shape[0]

    point_indices_to_connect = []
    iterated_points = []
    t1 = time.time()
    for i in range(full_list_length):
        if full_list[i] in iterated_points:
            continue
        iterated_points.append(full_list[i])
        for j in range(i+1, full_list_length):
            if numpy.linalg.norm(full_list[i] - full_list[j]) < eps:
                point_indices_to_connect.append((i,j))

    t2 = time.time()

    reduced_list = combine_points(full_list, point_indices_to_connect)
    t3 = time.time()

    t_connect = t2 - t1
    t_remove  = t3 - t2

    print('point connection time:', t_connect)
    print('point removal time:', t_remove)

    return full_list, point_indices_to_connect, reduced_list


def combine_points(numpy.ndarray full_list, list point_indices_to_connect):
    cdef numpy.ndarray reduced_list
    reduced_list = numpy.delete(
        full_list,
        [indices[1] for indices in point_indices_to_connect],
        axis=0,
    )

    return reduced_list
