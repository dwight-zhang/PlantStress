#---------------- import ------------------
from ngsolve import *
from netgen.csg import *
import netgen.gui
from netgen.meshing import ImportMesh
from netgen.meshing import FaceDescriptor
import numpy as np
from scipy.linalg import eigh

#---------------- meshing ------------------
ngmesh = ImportMesh('400.mesh')
mesh = Mesh(ngmesh)
Draw(mesh)
input("Mesh Done!")

#---------------- parameters ------------------
Pi = np.pi
nu = 0.5
Yg_pz = 0.3
Yg_cz = 1
Pout = 0.2
cf = CoefficientFunction(Yg_pz)
Yg = cf + (Yg_cz - cf) * IfPos(x**2/a**2 + y**2/b**2 + z**2/c**2 - 1, 0, 1)
Draw(Yg, mesh, "stiff")
input("Stiff Done!")

#---------------- FESpace ------------------      
orderset = 1
dimset = mesh.dim
fes = H1(mesh, order=orderset,   dim=dimset)

#-------------------- energy function --------------------------
def SolveLEPDE(fes, Yg, nu, Pout):
    # TeTfunctions
    u = fes.TrialFunction()
    v = fes.TestFunction()

    # bilinear form
    lame1 = Yg / (2 * (1 + nu))
    lame2 = Yg  * nu /( (1 + nu) * ( 1 - 2 * nu))
    u_grad = grad(u).Trace()
    v_grad = grad(v).Trace()
    u_div = Trace(u_grad * Id(dimset)).Compile()
    v_div = Trace(u_grad * Id(dimset)).Compile()
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += lame1 * InnerProduct(u_grad, v_grad) * ds
    a += lame2 * u_div * v_div * ds
    a.Assemble()

    # linear form
    nsurf = specialcf.normal(dimset)
    force = Pout * nsurf
    f = LinearForm(fes)
    f += force*v*ds
    f.Assemble()

    # initial gridfunction and solve equation
    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs())*f.vec

    return gfu

#------------------- deformation functions ---------------------
gfu = GridFunction(fes)
gfu = SolveLEPDE(fes, Yg, nu, Pout)

#------------------- strain functions ---------------------
u_grad = grad(gfu)
I = Id(dimset)
lame1 = Yg / (2 * (1 + nu))
lame2 = Yg  * nu /( (1 + nu) * ( 1 - 2 * nu))
# E = lame1 * (u_grad.trans + u_grad - 2 * I) + lame2 * InnerProduct((u_grad.trans - I),I ) * I
E = lame1 * (u_grad.trans + u_grad - 2 * I) 
gfw = GridFunction(fes)
for v in mesh.vertices:
    vp = np.array(v.point)
    p = mesh(vp[0],vp[1],vp[2])
    array_E = np.array(E(p)).T
    new_E = array_E.reshape(3,3)
    w = eigh(new_E, eigvals_only=True)
    for j in range(3):
        gfw.vec[v.nr][j] = w[j]
Draw(gfu, mesh, "deform")
Draw(gfw, mesh, "strain")
nsurf = specialcf.normal(dimset)
surf = (I - OuterProduct(nsurf, nsurf))* gfw
Draw(surf, mesh, "surf")
input("Strain done!")

