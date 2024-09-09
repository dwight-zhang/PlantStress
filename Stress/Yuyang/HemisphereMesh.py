# -*- coding: utf-8 -*-
from netgen.csg import *
from ngsolve import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri
from netgen.libngpy._meshing import ElementId2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

def plot_mesh_3d(mesh, filename):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # ��ȡ���ж��������
    vertices = np.array([[point.p[0], point.p[1], point.p[2]] for point in mesh.ngmesh.Points()])
    print(f"Vertices: {vertices}")

    # ��ȡ����Ԫ�صĶ�����������������������Ӧ��0��ʼ������ϵͳ
    elements = []
    for el in mesh.ngmesh.Elements2D():
        element_vertices = [v.nr - 1 for v in el.vertices]  # ��1����Ӧ��0��ʼ������
        elements.append(element_vertices)
    print(f"Elements: {elements}")

    # ȷ�������Ԫ�ز�Ϊ��
    if len(vertices) == 0 or len(elements) == 0:
        print("Error: Vertices or elements are empty.")
        return

    # ��������
    for element in elements:
        poly3d = [[vertices[vertex] for vertex in element]]
        ax.add_collection3d(Poly3DCollection(poly3d, facecolors='w', edgecolors='k', linewidths=1, alpha=0.5))

    ax.set_title("3D Mesh Visualization")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # ������ķ�Χ��ȷ���������
    max_range = np.array([vertices[:,0].max()-vertices[:,0].min(), 
                          vertices[:,1].max()-vertices[:,1].min(), 
                          vertices[:,2].max()-vertices[:,2].min()]).max() / 2.0

    mid_x = (vertices[:,0].max()+vertices[:,0].min()) * 0.5
    mid_y = (vertices[:,1].max()+vertices[:,1].min()) * 0.5
    mid_z = (vertices[:,2].max()+vertices[:,2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.savefig(filename)
    plt.close()

order = 3

geo = CSGeometry()
sph = Sphere(Pnt(0, 0, 0), 1).bc("sphere")
down = Plane(Pnt(0, 0, 0), Vec(0, 0, -1))
sphanddown = sph * down
geo.AddSurface(sph, sphanddown)
geo.NameEdge(sph, down, "down")

mesh = Mesh(geo.GenerateMesh(maxh=0.2))
mesh.Curve(order)

# ��ӡ���������Ϣ
print(f"Generated mesh: {mesh.ngmesh}")
print(f"Number of points: {len(list(mesh.ngmesh.Points()))}")
print(f"Number of 2D elements: {len(list(mesh.ngmesh.Elements2D()))}")

# ���û�����ά����ĺ���
plot_mesh_3d(mesh, "mesh_3d.png")

# ��ȡ Netgen �� Mesh ����
ngmesh = mesh.ngmesh

# ʹ�� Netgen �� Export ������������Ϊ .mesh �ļ�
ngmesh.Export("HemisphereMesh.mesh", "Neutral Format")

print("Mesh saved as HemisphereMesh.mesh")

# ʹ�� VTKOutput ������������Ϊ .vtu �ļ�
vtk = VTKOutput(mesh, filename="HemisphereMesh", subdivision=0)

    # ��ӽڵ�͵�Ԫ��Ϣ
vtk.Do(vb=BND)

print("Mesh saved as HemisphereMesh.vtu")
