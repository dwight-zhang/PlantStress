# geometric non-linear elasticity with Neo-Hooke hyperelastic material

#---------------- import netgen/ngsolve ------------------
import netgen.gui
from netgen.csg import *
from ngsolve import *

#--------------------- generate mesh ---------------------
geo = CSGeometry()
sphere = Sphere(Pnt(0,0,0),1)
plane  = Plane(Pnt(0,0,0), Vec(0,0,-1) ).bc("left")
geo.Add(sphere*plane.maxh(0.05))
ngmesh = geo.GenerateMesh()
mesh = Mesh(ngmesh)

#--------------------- parameters ------------------------
E, nu = 210, 0.2
mu  = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

#--------------------- FESpace ---------------------------
fes = H1(mesh, order=2, dirichlet="left", dim=mesh.dim)
# fes = VectorH1(mesh, order=2, dirichlet="left")

#-------------------- Neo-Hooke --------------------------
u  = fes.TrialFunction()
force = CoefficientFunction((0,0,1/2))
I = Id(mesh.dim)
F = I + Grad(u)
C = F.trans * F
E = 0.5 * (C-I)
def Pow(a, b):
    return a**b  # exp (log(a)*b)
def NeoHooke (C):
    return 0.5 * mu * (Trace(C-I) + 2*mu/lam * Pow(Det(C),-lam/2/mu) - 1)
factor = Parameter(1)
a = BilinearForm(fes, symmetric=False)
a += Variation(NeoHooke(C).Compile()*dx)
a += Variation((-factor * InnerProduct(force,u)).Compile()*dx)

#------------------- grid functions ---------------------
u = GridFunction(fes)
u.vec[:] = 0
res = u.vec.CreateVector()
w = u.vec.CreateVector()

#------------------- Newton iteration -------------------
for loadstep in range(10):
    # print loadstep number
    print ("loadstep", loadstep)
    factor.Set (loadstep+1)
    # Newton iteration energy 
    for it in range(5):
        print ("Newton iteration", it)
        print ("energy = ", a.Energy(u.vec))
        a.Apply(u.vec, res)
        a.AssembleLinearization(u.vec)
        inv = a.mat.Inverse(fes.FreeDofs())
        w.data = inv*res
        print ("err^2 = ", InnerProduct (w,res))
        u.vec.data -= w
    # draw dispalcement
    Draw(u, mesh, "displacement")
    SetVisualization(deformation=True)
    input("<press a key>")
