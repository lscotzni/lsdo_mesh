import lsdo_mesh.geometry_classes as lm
import numpy as np
import matplotlib.pyplot as plt

f               = lm.FFD()
a               = 5

ffd_point_1     = lm.Vertex(0.,0.)
ffd_point_2     = lm.Vertex(a,0.)
ffd_point_3     = lm.Vertex(a,a)
ffd_point_4     = lm.Vertex(0.,a)

ffd_points = [
    ffd_point_1,
    ffd_point_2,
    ffd_point_3,
    ffd_point_4
]

ffd_curve_1     = lm.Edge(ffd_point_1,ffd_point_2)
ffd_curve_2     = lm.Edge(ffd_point_2,ffd_point_3)
ffd_curve_3     = lm.Edge(ffd_point_3,ffd_point_4)
ffd_curve_4     = lm.Edge(ffd_point_4,ffd_point_1)

# ffd_face_1      = lm.Face([
#     ffd_point_1,
#     ffd_point_2,
#     ffd_point_3,
#     ffd_point_4,
# ], input_type='polygon')

ffd_face_1      = lm.Face([
    ffd_curve_1,
    ffd_curve_2,
    ffd_curve_3,
    ffd_curve_4],
    input_type='edges',
)

f.add_ffd_entity(ffd_face_1)

f.assemble()

print(vars(ffd_point_1)['properties']['x'])

control_points = np.zeros((2, 2, 2, 3))
# print(control_points)
# print(control_points.shape)

plt.figure(1)
for i in range(len(ffd_points)):
    if i == len(ffd_points) - 1:
        ps = list(vars(ffd_points[i])['properties'].values())
        pf = list(vars(ffd_points[0])['properties'].values())

        plt.plot([ps[0], pf[0]],[ps[1], pf[1]],'k')
    else:
        ps = list(vars(ffd_points[i])['properties'].values())
        pf = list(vars(ffd_points[i+1])['properties'].values())

        plt.plot([ps[0], pf[0]],[ps[1], pf[1]],'k')

# plt.show()