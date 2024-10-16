import netgen.gui
from ngsolve import *
from ngsolve.internal import visoptions
from netgen.meshing import ImportMesh

ngmesh = ImportMesh('cell_mesh_R1_opt.mesh')

from netgen.meshing import FaceDescriptor
isinplane = lambda el: all((ngmesh[v][2] < -27 for v in el.vertices))
fd_id = ngmesh.Add(FaceDescriptor(surfnr=2,domin=1,domout=0,bc=2))
for el in ngmesh.Elements2D():
    if isinplane(el):
        el.index = fd_id
ngmesh.SetBCName(fd_id-1, "new_bc")

mesh = Mesh(ngmesh)
Draw(mesh)

fes = H1(mesh, order = 1, dirichlet="new_bc")
u,v = fes.TnT()

a = BilinearForm(fes, symmetric=True)
a += grad(u).Trace() * grad(v).Trace()*ds
a.Assemble()

# external force
force = sin(x)
f = LinearForm(fes)
f += force*v*ds
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs())*f.vec
Draw(gfu, mesh, "u", autoscale=False)

gfu = GridFunction(fes)
gfu.Set(-sin(x), definedon=mesh.Boundaries("new_bc"))
r = f.vec.CreateVector()
r.data = f.vec - a.mat*gfu.vec
gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
Draw(gfu, mesh, "u", autoscale=False)