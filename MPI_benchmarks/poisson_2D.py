# Sovles -laplace(u)=1 on [0,1]^2
# with Dirichelet boundary u=0 
# exact solution 16*x*(1-x)*y*(1-y)
#
# Call with:
# mpirun -np 4 python ngs_poisson_2D.py
#
# Solver: 
# direct or CGSolver with preconditioner (Jacobi, BoomerAMG, BDDC/pardiso)

## imports
from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import unit_square

def master_print_line (a,b, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print ("\n"+a * 25,b.center(20),a * 25)
def master_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print (*args)
def worker_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank!=0:
        print (*args)
    
## initialize MPI
comm = MPI.COMM_WORLD
rank = comm.rank
np = comm.size
from ngsolve import __version__
master_print_line("=","LOAD OPTIONS")
master_print("NGSolve-"+__version__)
master_print(MPI.Get_library_version())
master_print("node =",MPI.Get_processor_name(),", num_process = ",np,", total_submesh = ", np-1)

## generate, distribute, and refine mesh & build H1-FESpace as usual
master_print_line("=","MESH & FEM SPACE")
t0 = MPI.Wtime()
if rank == 0:
    # master-proc generates mesh
    ngmesh = unit_square.GenerateMesh(maxh=0.1).Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
h=6 # refine-h
p=4 # order-p
for l in range(h):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
fes = H1(mesh, order=p, dirichlet=".*") # H1-P4
u,v = fes.TnT()
t1 = MPI.Wtime()
mesh_time = comm.allreduce(t1-t0, MPI.MAX)
tol_ele = comm.allreduce(mesh.GetNE(VOL), MPI.SUM) 
sum_local_dof = comm.allreduce(fes.ndof, MPI.SUM)
master_print("fes_space_name =",fes.name,", dim =",fes.dim)
master_print("mesh_reflevel =", h,", poly_order =",p)
worker_print("submesh", rank, ", local_ele =", mesh.GetNE(VOL),", local_dof =", fes.ndof)
master_print("total_eles =", tol_ele, ", global_dofs =", fes.ndofglobal) 
master_print("sum_local_dofs =", sum_local_dof, ", avg_local_dofs = ", int(sum_local_dof/(np-1)))
master_print("[max_wall_time =",mesh_time,"s]")


## BLF and preconditoners
master_print_line("=","MAT/RHS ASSEMBLE")
t0 = MPI.Wtime()
a = BilinearForm(grad(u)*grad(v)*dx)
# pre = Preconditioner(a, type="direct", inverse="masterinverse") # direct solve with mumps
# pre = Preconditioner(a, type="direct", inverse="sparsecholesky") # direct solve with sparsecholesky
# pre = Preconditioner(a, type="direct", inverse="pardiso") # direct solve with pardiso
# pre = Preconditioner(a, type="local") # Jacobi 
# pre = Preconditioner(a, type="hypre") # BoomerAMG (use only for order 1)
pre = Preconditioner(a, type="bddc", inverse="pardiso") # BDDC + pardiso for coarse matrix
master_print("pre =",pre.name)
a.Assemble()
pre.Update()
# RHS does not change either!
f = LinearForm(32 * (y*(1-y)+x*(1-x))*v*dx).Assemble()
t1 = MPI.Wtime()
assemble_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("[max_wall_time =",assemble_time,"s]")

## solve the equation
master_print_line("=","CGSolver")

from ngsolve.krylovspace import CGSolver
gfu = GridFunction(fes)
inv = CGSolver(a.mat, pre.mat, printing=False, maxiter=5000, tol=1e-8) # use CGsolver with precondition pre
gfu.vec.data = inv*f.vec
t1 = MPI.Wtime()
solve_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("absolute_error = ", inv.errors[len(inv.errors)-1])
master_print("[max_wall_time =",solve_time,"s]")

# check solution 
master_print_line("=","CHECK SOLUTION")
t0 = MPI.Wtime()
exact = 16*x*(1-x)*y*(1-y)
error = Integrate ( (gfu-exact)*(gfu-exact) , mesh)
t1 = MPI.Wtime()
check_time = comm.allreduce(t1-t0, MPI.MAX)
master_print ("L2-error =", error)
master_print("[max_wall_time =",check_time,"s]")

# visual the sub-domains obtained by the automatic partitioning
fesL2 = L2(mesh, order=0)
gfrank = GridFunction(fesL2)
gfrank.vec.local_vec[:] = comm.rank

# output by pickle
import os
import pickle
output_time = 0
do_pickle = False
if do_pickle:
    output_path = os.getcwd()+"../pickles/poisson_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT PICKLE")
    t0 = MPI.Wtime()
    netgen.meshing.SetParallelPickling(True)
    pickle.dump((gfrank, gfu),open(output_path+"/pickleout"+str(comm.rank), "wb"))
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",output_time,"s]")
    
# output by vtk
do_vtk = False
if do_vtk:
    output_path = os.getcwd()+"../vtks/poisson_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT VTK")
    t0 = MPI.Wtime()
    vtk = VTKOutput(mesh, coefs=[gfu,gfrank], names=["u","rank"], filename=output_path+"/vtkout", subdivision=2)
    vtk.Do()
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",vtk_time,"s]")

# wall time 
wall_time = mesh_time + assemble_time + solve_time + check_time + output_time
master_print_line("=","TOTAL WALL TIME")
master_print("[tol_wall_time =",wall_time,"s]")