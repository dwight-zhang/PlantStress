#---------------- import ------------------
from ngsolve import *
from netgen.csg import *
import netgen.gui
from netgen.meshing import ImportMesh
from netgen.meshing import FaceDescriptor
import numpy as np
from scipy.integrate import quad

#---------------- meshing ------------------
geo = CSGeometry()
sph = Sphere(Pnt(0, 0, 0.2), 1).bc("sph")
down = Plane(Pnt(0, 0, 0), Vec(0, 0, -1)).bc("down")
top = Plane(Pnt(0, 0, 1.1), Vec(0, 0, -1)).bc("top")
sphtop = sph * top
sphdown = sph * down
pt = sphdown - sphtop 
geo.AddSurface(sph, pt)
geo.AddSurface(down, pt)
geo.AddSurface(top, pt)
ngmesh = geo.GenerateMesh(maxh=0.08)
mesh = Mesh(ngmesh)
mesh.Curve(2)
Draw(mesh)

#---------------- parameters ------------------
Pi = np.pi
orderset = 2
dimset = mesh.dim
v = 0.5
Yg = 1
Pout = 1

#---------------- FESpace ------------------      
fes = H1(mesh, order=orderset,  dirichlet="down|top", dim=dimset)
u = fes.TrialFunction()
nsurf = -specialcf.normal(3)

#-------------------- energy function --------------------------
a = BilinearForm(fes, symmetric=False, check_unused=False)
u_grad = grad(u).Trace()
Sgradu = u_grad - u_grad * OuterProduct(nsurf, nsurf)
Eu = 1/2*(Sgradu.trans * Sgradu - Id(dimset))
Cg = CoefficientFunction(((Yg,v*Yg,0),(Yg,v*Yg,0),(0,0,Yg*(v-1)/2)),dims=(3,3))
Phi = Eu.trans * (Cg * Eu)
a += Variation(Trace(Phi).Compile()*ds)
ngradu = u_grad.trans * nsurf
beta = Pout*nsurf
a += Variation((ngradu-beta)*(ngradu-beta)*ds(definedon=fes.mesh.Boundaries("sph")))

#------------------- grid functions ---------------------
u = GridFunction(fes)
du = GridFunction(fes)
res = u.vec.CreateVector()
# ini = CoefficientFunction([sin(x*x*y*y)*z*(1.1-z)*nsurf if bc == 'sph' else (0,0,0) for bc in mesh.GetBoundaries()])
ini = CoefficientFunction([sin(x*x*y*y)*z*(1.1-z)*nsurf])
u.Set(ini,definedon=fes.mesh.Boundaries("sph"))
Draw(u, mesh, 'displacement')
input("initial value setted!")

#------------------- damped Newton ---------------------
step = 0.2
for i in range(50):
   print ("Damped Newton iteration", i+1)
   print ("  energy = ", a.Energy(u.vec))
   a.Apply(u.vec, res) 
   a.AssembleLinearization(u.vec)
   du.vec.data = a.mat.Inverse(fes.FreeDofs()) * res
   r =  InnerProduct(du.vec, res)
   print ("  err^2 =", abs(r))
   u.vec.data -= step * du.vec.data
   Draw(u, mesh, 'displacement')
   input("  iteration done!")
   if abs(r) < 1e-6:
       break

