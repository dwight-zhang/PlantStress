# -*- coding: utf-8 -*-

#---------------- import packages ------------------

from ngsolve import *
from netgen.csg import *
from netgen.meshing import ImportMesh
from netgen.meshing import FaceDescriptor
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.integrate import quad

#---------------- draw mesh ------------------

# def plot_mesh_3d(mesh, filename):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     vertices = np.array([[point.p[0], point.p[1], point.p[2]] for point in mesh.ngmesh.Points()])

#     elements = []
#     for el in mesh.ngmesh.Elements2D():
#         element_vertices = [v.nr - 1 for v in el.vertices]
#         elements.append(element_vertices)

#     if len(vertices) == 0 or len(elements) == 0:
#         print("Error: Vertices or elements are empty.")
#         return

#     for element in elements:
#         poly3d = [[vertices[vertex] for vertex in element]]
#         ax.add_collection3d(Poly3DCollection(poly3d, facecolors='w', edgecolors='k', linewidths=1, alpha=0.5))

#     ax.set_title("3D Mesh Visualization")
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')

#     max_range = np.array([vertices[:,0].max()-vertices[:,0].min(), 
#                           vertices[:,1].max()-vertices[:,1].min(), 
#                           vertices[:,2].max()-vertices[:,2].min()]).max() / 2.0

#     mid_x = (vertices[:,0].max()+vertices[:,0].min()) * 0.5
#     mid_y = (vertices[:,1].max()+vertices[:,1].min()) * 0.5
#     mid_z = (vertices[:,2].max()+vertices[:,2].min()) * 0.5

#     ax.set_xlim(mid_x - max_range, mid_x + max_range)
#     ax.set_ylim(mid_y - max_range, mid_y + max_range)
#     ax.set_zlim(mid_z - max_range, mid_z + max_range)

#     plt.savefig(filename)
#     plt.close()



#--------------------- load mesh ---------------------

# 导入网格
ngmesh = ImportMesh('HemisphereMesh.mesh')

#最低的一层定义为固定边界
isinplane = lambda el: any((ngmesh[v][2] <= 0 for v in el.vertices))
fd_id = ngmesh.Add(FaceDescriptor(surfnr=2,domin=1,domout=0,bc=2))

for el in ngmesh.Elements2D():
  if isinplane(el):
    el.index = fd_id
  else:
    el.index = fd_id - 1

ngmesh.SetBCName(fd_id, "down")

mesh = Mesh(ngmesh)

# plot_mesh_3d(mesh, "mesh_3d.png")

#---------------- parameters ------------------

orderset = 2

dimset = mesh.dim

v = 0.5
Yg = 1
Yf = 1
Pi = np.pi

Pout = 1

alpha0 = 100

def rho(theta):
  return np.sin(theta) + 10000 * np.cos(theta)

def integrand0(theta):
    return rho(theta)

# 计算实部和虚部的 Fourier 系数
rho0value, _ = quad(lambda theta: integrand0(theta).real, 0, np.pi)
rho0value = rho0value*2

# 定义被积函数
def integrand1(theta):
    return rho(theta) * np.exp(-2j * theta)

# 计算实部和虚部的 Fourier 系数
rho1value, _ = quad(lambda theta: integrand1(theta).real, 0, np.pi)
rho1bvalue, _ = quad(lambda theta: integrand1(theta).imag, 0, np.pi)
rho1value = rho1value*2
rho1bvalue = - rho1bvalue*2

# 定义被积函数
def integrand2(theta):
    return rho(theta) * np.exp(-4j * theta)

# 计算实部和虚部的 Fourier 系数
rho2value, _ = quad(lambda theta: integrand2(theta).real, 0, np.pi)
rho2bvalue, _ = quad(lambda theta: integrand2(theta).imag, 0, np.pi)
rho2value = rho2value*2
rho2bvalue = - rho2bvalue*2

#--------------------- FESpace ---------------------------

fes = H1(mesh, order=orderset, dirichlet="down", dim=dimset)

#-------------------- definite operators --------------------------

u = fes.TrialFunction()

nsurf = specialcf.normal(3) # 定义求法向量操作

#-------------------- external parameter functions --------------------------

P = Parameter(1)
P.Set(Pout) # 将膨压力设置为均匀的定值

rho0 = Parameter(2)
rho0.Set(rho0value)

rho1 = Parameter(3)
rho1.Set(rho1value)

rho1b = Parameter(4)
rho1b.Set(rho1bvalue)

rho2 = Parameter(5)
rho2.Set(rho2value)

rho2b = Parameter(6)
rho2b.Set(rho2bvalue)

#-------------------- energy function --------------------------

u_grad = grad(u).Trace()

Sgradu = InnerProduct(nsurf, nsurf)*u_grad - u_grad*OuterProduct(nsurf, nsurf)

Eu = 1/2*(Sgradu.trans*Sgradu - Id(dimset))

code='''
for i in range({0}):
  for j in range({0}):
    tempmatrix = [[0 for j in range({0})] for i in range({0})]
    for i1 in range({0}):
      for j1 in range({0}):
        if i1 == i and j1 == j:
          tempmatrix[i1][j1] = 1
        else:
          tempmatrix[i1][j1] = 0
    codein1 = "E{{0}}{{1}}=CoefficientFunction(tempmatrix)"
    exec(codein1.format(i,j))
    codein2 = "global E{{0}}{{1}}"
    exec(codein2.format(i,j))
'''
exec(code.format(dimset))

def Cmatrix(rho0, rho1,rho1b,rho2,rho2b):
  Cmat = Yg*(1*(globals()['E00']+globals()['E11']) + v*(globals()['E10']+globals()['E01']) + (1-v)/2*globals()['E22']) + Pi*Yf/16*((3*rho0+rho2+4*rho1)*globals()['E00']\
  +(3*rho0+rho2-4*rho1)*globals()['E11']+(rho0-rho2)*(globals()['E01']+globals()['E10']+globals()['E22'])+(2*rho1b+rho2b)*(globals()['E02']+globals()['E20'])\
  +(2*rho1b-rho2b)*(globals()['E12']+globals()['E21']))
  return Cmat

Phi = Eu.trans * (Cmatrix(rho0, rho1, rho1b, rho2, rho2b) * Eu)
  
a = BilinearForm(fes, symmetric=False, check_unused=False)
a += Variation(Trace(Phi).Compile()*ds)

# 定义外力能量
# Pn = P/3 * InnerProduct(nsurf, u)
# a += Variation(Pn.Compile()*ds)

B = u_grad * nsurf - P/3 * nsurf
BF = alpha0 * InnerProduct(B, B)
a += Variation(BF.Compile()*ds)

#------------------- grid functions ---------------------

# 定义牛顿法求解过程
u = GridFunction(fes)
du = u.vec.CreateVector()
res = u.vec.CreateVector()
u.Set(CoefficientFunction((x, 1, 1)))

for i in range(10):
    print ("Newton iteration", i+1)
    print ("energy = ", a.Energy(u.vec))
    a.Apply(u.vec, res)  
    a.AssembleLinearization(u.vec)
    du.data = a.mat.Inverse(fes.FreeDofs()) * res
    print ("err^2=", InnerProduct(du, res))
    u.vec.data -= du.data

# VTK输出
print("Starting VTK output.")
vtk = VTKOutput(mesh, coefs=[u], names=["uoutput"], filename="result", subdivision=1)

# 添加节点和单元信息
vtk.Do(vb=BND)
print("VTK output completed. Check 'result.vtu' file.")
