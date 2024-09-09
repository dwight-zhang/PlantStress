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

    # 提取所有顶点的坐标
    vertices = np.array([[point.p[0], point.p[1], point.p[2]] for point in mesh.ngmesh.Points()])
    print(f"Vertices: {vertices}")

    # 提取所有元素的顶点索引，并调整索引以适应从0开始的索引系统
    elements = []
    for el in mesh.ngmesh.Elements2D():
        element_vertices = [v.nr - 1 for v in el.vertices]  # 减1以适应从0开始的索引
        elements.append(element_vertices)
    print(f"Elements: {elements}")

    # 确保顶点和元素不为空
    if len(vertices) == 0 or len(elements) == 0:
        print("Error: Vertices or elements are empty.")
        return

    # 绘制网格
    for element in elements:
        poly3d = [[vertices[vertex] for vertex in element]]
        ax.add_collection3d(Poly3DCollection(poly3d, facecolors='w', edgecolors='k', linewidths=1, alpha=0.5))

    ax.set_title("3D Mesh Visualization")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # 计算轴的范围以确保网格居中
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

# 打印网格对象信息
print(f"Generated mesh: {mesh.ngmesh}")
print(f"Number of points: {len(list(mesh.ngmesh.Points()))}")
print(f"Number of 2D elements: {len(list(mesh.ngmesh.Elements2D()))}")

# 调用绘制三维网格的函数
plot_mesh_3d(mesh, "mesh_3d.png")

# 获取 Netgen 的 Mesh 对象
ngmesh = mesh.ngmesh

# 使用 Netgen 的 Export 方法保存网格为 .mesh 文件
ngmesh.Export("HemisphereMesh.mesh", "Neutral Format")

print("Mesh saved as HemisphereMesh.mesh")

# 使用 VTKOutput 方法保存网格为 .vtu 文件
vtk = VTKOutput(mesh, filename="HemisphereMesh", subdivision=0)

    # 添加节点和单元信息
vtk.Do(vb=BND)

print("Mesh saved as HemisphereMesh.vtu")
