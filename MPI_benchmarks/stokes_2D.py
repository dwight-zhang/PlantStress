# Solves - laplace(u) - grad(p) = f
#        div(u)  = 0
# with inlet boundary u - 1.5*4*y*(0.4-y)/(-.41*0.41)
#
# call with:
# mpirun -np 4 python ngs_stokes_2D.py
#
# Solver:
# Diagnal block matrix preconditioner
#
# FEM:
# Taylor-Hood elements

from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import SplineGeometry
from ngsolve.krylovspace import CGSolver
from ngsolve.krylovspace import MinRes

def master_print_line (a,b, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print ("\n"+a * 25,b.center(20),a * 25)
def master_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print (*args)
def worker_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank!=0:
        print (*args)

h = 0 # mesh size (<=3)
k = 2 # poly order (k>=2)
# pre_use = 'diag-exact'
pre_use = 'AL-exact'
# pre_use = 'diag-inexact'
# pre_use = 'AL-inexact'

## initialize MPI
comm = MPI.COMM_WORLD
rank = comm.rank
np = comm.size
from ngsolve import __version__
master_print_line("=","LOAD OPTIONS")
master_print("NGSolve-"+__version__)
master_print(MPI.Get_library_version())
print("node =",MPI.Get_processor_name(),", num_process = ",np,", total_submesh = ", np-1)

# mesh & Taylor-hood elements (P^k)^2*P^(k-1)
master_print_line("=","MESH & FEM SPACE")
t0 = MPI.Wtime()
geo = SplineGeometry()
geo.AddRectangle ((0, 0), (2, 0.41), bcs = ("wall", "outlet", "wall", "inlet"))
geo.AddCircle ((0.2, 0.2), r=0.05, leftdomain=0, rightdomain=1, bc="cyl")
if rank ==0:
    ngmesh = geo.GenerateMesh(maxh=0.02)
    ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
    ngmesh.SetGeometry(geo)
for l in range(h):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
V = VectorH1(mesh, order=k, dirichlet="wall|inlet|cyl")
Q = H1(mesh, order=k-1)
X = FESpace([V,Q])
t1 = MPI.Wtime()
mesh_time = comm.allreduce(t1-t0, MPI.MAX)
tol_ele = comm.allreduce(mesh.GetNE(VOL), MPI.SUM) 
sum_local_dof = comm.allreduce(X.ndof, MPI.SUM)
master_print("mesh_refine_level =", h,", poly_order =",k)
worker_print("submesh", rank, ", local_ele =", mesh.GetNE(VOL),", local_dof =", X.ndof)
master_print("total_eles =", tol_ele, ", global_dofs =", X.ndofglobal) 
master_print("sum_local_dofs =", sum_local_dof, ", avg_local_dofs = ", int(sum_local_dof/(np-1)))
master_print("[max_wall_time =",mesh_time,"s]")


## BLF and preconditoners
master_print_line("=","MAT/RHS ASSEMBLE")
t0 = MPI.Wtime()
u,v = V.TnT()
p,q = Q.TnT()
a = BilinearForm(V)
a += InnerProduct(Grad(u),Grad(v))*dx
b = BilinearForm(trialspace=V, testspace=Q)
b += div(u)*q*dx
m = BilinearForm(Q)
m += p*q*dx
f = LinearForm(V)
f += CoefficientFunction((0,x-0.5)) * v * dx
g = LinearForm(Q)
# auxiliary preconditioner
if pre_use == 'diag-inexact':
    ap = Preconditioner(a, type='bddc', inverse='pardiso')
if pre_use == 'diag-inexact' or pre_use == 'AL-inexact':
    mp = Preconditioner(m, type='bddc', inverse='pardiso')
f.Assemble()
g.Assemble()
m.Assemble()
a.Assemble()
b.Assemble()
K = BlockMatrix ([[a.mat, b.mat.T], [b.mat, None]])
## block diagonal preconditioning matrices (diag-pre)
if pre_use == 'diag-exact':
    minv = m.mat.Inverse(inverse='sparsecholesky')
    ainv = a.mat.Inverse(V.FreeDofs(), inverse='sparsecholesky')
    C = BlockMatrix ([[ainv, None], [None, minv]])
elif pre_use == 'diag-inexact':
    C = BlockMatrix ([[ap.mat, None], [None, mp.mat]])
## the Augmented Lagrangian method as a preconditioner (AL-pre)
elif pre_use == 'AL-exact':
    eps = 1e-9
    al = BilinearForm(V)
    al += (InnerProduct(Grad(u),Grad(v))+eps*div(u)*div(v))*dx
    al.Assemble()
    minv = m.mat.Inverse(inverse='sparsecholesky')
    alinv = al.mat.Inverse(V.FreeDofs(), inverse='sparsecholesky')
    C = BlockMatrix ([[alinv, None], [- eps * b.mat @ alinv,eps * minv]])
elif pre_use == 'AL-inexact':
    eps = 0.5
    al = BilinearForm(V)
    al += (InnerProduct(Grad(u),Grad(v))+eps*div(u)*div(v))*dx
    alp = Preconditioner(al,  type='bddc', inverse='pardiso')
    al.Assemble()
    C = BlockMatrix ([[alp.mat, None], [- eps * b.mat @ alp.mat,eps * mp.mat]])
t1 = MPI.Wtime()
assemble_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("pre =",pre_use)
master_print("[max_wall_time =",assemble_time,"s]")

# Set boundary & MinRes solve
master_print_line("=","BDRY & MinRes")
t0 = MPI.Wtime()
gfu = GridFunction(V, name="u")
gfp = GridFunction(Q, name="p")
uin = CoefficientFunction( (1.5*4*y*(0.41-y)/(0.41*0.41), 0) )
gfu.Set(uin, definedon=mesh.Boundaries("inlet"))
rhs = BlockVector ([f.vec, g.vec] )
sol = BlockVector ([gfu.vec, gfp.vec])
MinRes(mat=K, pre=C, rhs=rhs, sol=sol, initialize=False, maxsteps=50000, printrates=comm.rank==0, tol=1e-08)
t1 = MPI.Wtime()
solve_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("[max_wall_time =",solve_time,"s]")

# visual the sub-domains obtained by the automatic partitioning
fesL2 = L2(mesh, order=0)
gfrank = GridFunction(fesL2)
gfrank.vec.local_vec[:] = comm.rank

# output by pickle
output_time = 0
import os
import pickle
do_pickle = False
if do_pickle:
    output_path = os.getcwd()+"/../pickles/stokes_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT PICKLE")
    t0 = MPI.Wtime()
    netgen.meshing.SetParallelPickling(True)
    pickle.dump((gfrank, gfu, gfp),open(output_path+"/pickleout"+str(comm.rank), "wb"))
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",output_time,"s]")
    
# output by vtk
do_vtk = False
if do_vtk:
    output_path = os.getcwd()+"/../vtks/stokes_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT VTK")
    t0 = MPI.Wtime()
    vtk = VTKOutput(mesh, coefs=[gfrank,gfu,gfp], names=["rank","u","p"], filename=output_path+"/vtkout", subdivision=2)
    vtk.Do()
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",vtk_time,"s]")
    
# wall time 
wall_time = mesh_time + assemble_time + solve_time + output_time
master_print_line("=","TOTAL WALL TIME")
master_print("[tol_wall_time =",wall_time,"s]")