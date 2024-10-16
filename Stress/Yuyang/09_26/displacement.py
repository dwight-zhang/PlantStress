#---------------- import ------------------
from ngsolve import *
from netgen.csg import *
import netgen.gui
from netgen.meshing import ImportMesh
from netgen.meshing import FaceDescriptor
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

#---------------- meshing ------------------
R = 1
pz_down = 0.0
cz_down = 0.8
geo = CSGeometry()
sph1 = Sphere(Pnt(0, 0, 0), R).bc("pz")
sph2 = Sphere(Pnt(0, 0, 0), R).bc("cz")
down = Plane(Pnt(0, 0, pz_down), Vec(0, 0, -1)).bc("down")
top = Plane(Pnt(0, 0, cz_down), Vec(0, 0, -1)).bc("top")
sph1down = sph1 * down
sph2top = sph2 * top
pt = sph1down - sph2top
geo.AddSurface(sph1, pt)
geo.AddSurface(down, pt)
geo.AddSurface(sph2, sph2top)
ngmesh = geo.GenerateMesh(maxh=0.1)
mesh = Mesh(ngmesh)
mesh.Curve(2)
Draw(mesh)
input("Mesh Done!")

#---------------- parameters ------------------
Pi = np.pi
nu = 0.5
pzofcz = 0.25
Yg_cz = 1
Pout = 0.2
Yg = CoefficientFunction([ pzofcz * Yg_cz if bc == "pz" else Yg_cz for bc in mesh.GetBoundaries()])
Draw(Yg, mesh, "stiff")
input("Stiff Done!")

#---------------- FESpace ------------------      
orderset = 2
dimset = mesh.dim
fes = H1(mesh, order=orderset,  dirichlet="down", dim=dimset)

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

