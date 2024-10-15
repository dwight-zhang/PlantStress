#---------------- import ------------------
from ngsolve import *
from netgen.csg import *
import netgen.gui
from netgen.meshing import ImportMesh
from netgen.meshing import FaceDescriptor
import numpy as np
from scipy.integrate import quad

#---------------- meshing ------------------
R = 1
pz_down = 0.0
cz_down = 0.75
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
ngmesh = geo.GenerateMesh(maxh=0.05)
mesh = Mesh(ngmesh)
mesh.Curve(2)
Draw(mesh)
input("Mesh Done!")

#---------------- parameters ------------------
Pi = np.pi
v = 0.5
pzofcz = 0.3
Yg_cz = 1e-1
Pout = 0.1
Yg = CoefficientFunction([ pzofcz * Yg_cz if bc == "pz" else Yg_cz for bc in mesh.GetBoundaries()])

#---------------- FESpace ------------------      
orderset = 2
dimset = mesh.dim
fes = H1(mesh, order=orderset,  dirichlet="down", dim=dimset)
u = fes.TrialFunction()
v = fes.TestFunction()
nsurf = specialcf.normal(dimset)

#--------------- Bilinear and liner form ----------------
# penalty = 1e4
u_grad = grad(u).Trace()
Eu = 1/2*(u_grad.trans +  u_grad) - Id(dimset)
Cg = CoefficientFunction(((Yg,v*Yg,0),
                          (v*Yg,Yg,0),
                          (0,0,Yg*(1-v)/2)),
                          dims = (dimset,dimset))
Phi = Eu.trans * (Cg * Eu)
ngradu = u_grad * nsurf
beta = Pout * nsurf
Force = CoefficientFunction((0,0,1))
a = BilinearForm(fes, symmetric=False, check_unused=False)
a += Variation(Trace(Phi).Compile()*ds(definedon=fes.mesh.Boundaries("pz|cz")))
a += Variation((InnerProduct(- Pout * nsurf,u))*ds(definedon=fes.mesh.Boundaries("pz|cz")))

#------------------- grid functions ---------------------
u = GridFunction(fes)
du = GridFunction(fes)
res = u.vec.CreateVector()
ini = CoefficientFunction(Pout*(z-pz_down)*nsurf)
u.Set(ini,definedon=fes.mesh.Boundaries("pz|cz"))
Draw(u, mesh, 'displacement')
input("Initiate Done!")

#------------------- FEM Solver ---------------------
step = 1
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


u_grad = grad(u).Trace()
Eu = 1/2*(u_grad.trans * u_grad - Id(dimset))
print(u_grad)
Draw( Eu, mesh, 'strain')
input("Done!")
