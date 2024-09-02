# Solves lid-driven cavity problem:
#     - nu * laplace(u) + u * grad(u) + grad(p) = f
#                                     - div(u)  = 0
# with dirichlet boundary u = ( ,0) on top,
#      domain = unit_square
#
# call with:
# mpirun -np 4 python steady_navierstokes_2D.py
#
# Solver:
# Diagnal block matrix preconditioner
#
# FEM:
# Taylor-Hood elements

h = 1                       # refine-h  
k = 6                       # (P^k)^2*P^(k-1)  >= 2
Re = 1000                   # Reynold number 
picard_steps = 5            # picard steps
maxits = 200                # picard+newton steps 
dampfactor = 1              # damped nonlinear iteration (0< value <=1)
nonlinear_tol = 1e-8        # nonlinear iteration steps
gmr_tol = 1e-8              # GMRES - tolerance (maybe sensitive)
# pre_use = 'none'            # direct solve
pre_use = 'al-exact'        # GMRES - pre use AL-exact
# pre_use = 'diag-exact'       # GMRES - pre use diag-exact
# pre_use = 'diag-inexact'    # need debug!
if pre_use == 'none':
    gam = 0.0 
else:
    gam = 1.0

## imports
from mpi4py import MPI
from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.krylovspace import CGSolver
from ngsolve.krylovspace import GMRes
from toolbox import LinMat,SymmetricGS

def master_print_line (a,b, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print ("\n"+a * 25,b.center(20),a * 25)
def master_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank==0:
        print (*args)
def worker_print (*args, comm=MPI.COMM_WORLD):
    if comm.rank!=0:
        print (*args)
## parameters
nu = Parameter(1/Re)
master_print("Re =", Re,", nu=",nu)
master_print("mesh_refine_level =", h,", poly_order =",k)

## initialize MPI
comm = MPI.COMM_WORLD
rank = comm.rank
np = comm.size
from ngsolve import __version__
master_print_line("=","LOAD OPTIONS")
master_print("NGSolve-"+__version__)
master_print(MPI.Get_library_version())
master_print("node =",MPI.Get_processor_name(),
             ", num_process = ",np,
             ", total_submesh = ", np-1)

## generate, distribute, and refine mesh & build H1-FESpace as usual
master_print_line("=","MESH & FEM SPACE")
t0 = MPI.Wtime()
if rank == 0:
    # master-proc generates mesh
    ngmesh = unit_square.GenerateMesh(maxh=0.05).Distribute(comm) 
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
for l in range(h):
    ngmesh.Refine()
mesh = Mesh(ngmesh)
V = VectorH1(mesh, order=k, dirichlet="bottom|right|top|left")
Q = H1(mesh, order=k-1)
# N = NumberSpace(mesh)
X = FESpace([V,Q])
t1 = MPI.Wtime()
mesh_time = comm.allreduce(t1-t0, MPI.MAX)
tol_ele = comm.allreduce(mesh.GetNE(VOL), MPI.SUM) 
sum_local_dof = comm.allreduce(X.ndof, MPI.SUM)
worker_print("submesh", rank, ", local_ele =", mesh.GetNE(VOL),
             ", local_dof =", X.ndof)
master_print("total_eles =", tol_ele, ", global_dofs =", X.ndofglobal) 
master_print("sum_local_dofs =", sum_local_dof, 
             ", avg_local_dofs = ", int(sum_local_dof/(np-1)))
master_print("[max_wall_time =",mesh_time,"s]")


# newton linearization iteration
master_print_line("=","MAT/RHS ASSEMBLE")
t0 = MPI.Wtime()
# boundary conditions
gfu = GridFunction(X)
gfu.components[0].Set(CoefficientFunction((1,0)),
                      definedon=mesh.Boundaries("top"))
# BLF initial
(u,p),(v,q) = X.TnT()
ns = BilinearForm(X)
ns += (nu*InnerProduct(grad(u),grad(v))
           +InnerProduct(grad(u)*u,v)
           +gam * div(u) * div(v)
           -div(u)*q-div(v)*p)*dx
picard = BilinearForm(X)
picard += (nu*InnerProduct(grad(u),grad(v))
           +InnerProduct(grad(u)*gfu.components[0],v)
           +gam * div(u) * div(v)
           -div(u)*q-div(v)*p)*dx
u,v = V.TnT()
p,q = Q.TnT()
picard_a = BilinearForm(V)
picard_a += (nu*InnerProduct(grad(u),grad(v))
           +InnerProduct(grad(u)*gfu.components[0],v)
           +gam * div(u) * div(v))*dx
newton_a = BilinearForm(V)
newton_a += (nu*InnerProduct(grad(u),grad(v))
           +InnerProduct(grad(u)*gfu.components[0],v)
           +InnerProduct(grad(gfu.components[0])*u,v)
           +gam * div(u) * div(v))*dx
b = BilinearForm(trialspace=V, testspace=Q)
b += -div(u)*q*dx
m = BilinearForm(Q)
m += p*q*dx
if pre_use == 'al-exact' or pre_use == 'diag-exact':
    b.Assemble(reallocate=True)
    m.Assemble(reallocate=True)
    minv = m.mat.Inverse(inverse='pardiso')
if pre_use == 'diag-inexact':
    minv = Preconditioner(m, type='local')
    b.Assemble(reallocate=True)
    m.Assemble()
t1 = MPI.Wtime()
assemble_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("pre =",pre_use)
master_print("[max_wall_time =",assemble_time,"s]")

# Simple Newton iteration 
master_print_line("=","Picard"+str(picard_steps)+"+Newton")
master_print("dampfactor =",dampfactor,
             ", picard_steps =", picard_steps)
res = gfu.vec.CreateVector()
du = gfu.vec.CreateVector()
for it in range(maxits):
    master_print_line('-',"Iteration "+str(it))
    res[:] = 0.0
    ns.Apply(gfu.vec,res)
    if it < picard_steps:
        picard.Assemble(reallocate=True)
        if pre_use == 'al-exact':
            picard_a.Assemble(reallocate=True)
            picard_a.mat.Update()
            picard_ainv = picard_a.mat.Inverse(V.FreeDofs(), 
                                               inverse='pardiso')
            upblock = (1/Re+gam) * picard_ainv @ b.mat.T @ minv
            AL = BlockMatrix ([[picard_ainv, upblock], 
                               [None, -(1/Re+gam) * minv]])
            AL_lin = LinMat(AL, [picard_a.mat.col_pardofs,
                                 b.mat.col_pardofs]) 
            GMRes(picard.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'diag-exact':
            picard_a.Assemble(reallocate=True)
            picard_a.mat.Update()
            picard_ainv = picard_a.mat.Inverse(V.FreeDofs(), 
                                               inverse='pardiso')
            AL = BlockMatrix ([[picard_ainv, None], 
                               [None, -(1/Re+gam) * minv]])
            AL_lin = LinMat(AL, [picard_a.mat.col_pardofs,
                                 b.mat.col_pardofs]) 
            GMRes(picard.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'diag-inexact':
            picard_ainv = Preconditioner(picard_a, type='local')
            picard_a.Assemble(reallocate=True)
            picard_a.mat.Update()
            AL = BlockMatrix ([[picard_ainv, None], 
                               [None, -(1/Re+gam) * minv]])
            AL_lin = LinMat(AL, [picard_a.mat.col_pardofs,
                                 b.mat.col_pardofs]) 
            GMRes(picard.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'none':
            du.data = picard.mat.Inverse(X.FreeDofs(),
                                         inverse='pardiso') * res
        # correction for gfu
        gfu.vec.data -= du
    else:
        ns.AssembleLinearization(gfu.vec,reallocate=True)
        if pre_use == 'al-exact':
            newton_a.Assemble(reallocate=True)
            newton_a.mat.Update()
            newton_ainv = newton_a.mat.Inverse(V.FreeDofs(), 
                                               inverse='pardiso')
            upblock = (1/Re+gam) * newton_ainv @ b.mat.T @ minv
            AL = BlockMatrix ([[newton_ainv, upblock], 
                               [None, -(1/Re+gam) * minv]])
            AL_lin = LinMat(AL, [newton_a.mat.col_pardofs,
                                 b.mat.col_pardofs]) 
            GMRes(ns.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'diag-exact':
            newton_a.Assemble(reallocate=True)
            newton_a.mat.Update()
            newton_ainv = newton_a.mat.Inverse(V.FreeDofs(), 
                                               inverse='pardiso')
            AL = BlockMatrix ([[newton_ainv, None], 
                               [None, -(1/Re+gam) * minv]])
            AL_lin = LinMat(AL, [newton_a.mat.col_pardofs, 
                                 b.mat.col_pardofs]) 
            GMRes(ns.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'diag-inexact':
            picard_ainv = Preconditioner(newton_a, type='local')
            newton_a.Assemble(reallocate=True)
            newton_a.mat.Update()
            AL = BlockMatrix ([[newton_ainv, None], 
                               [None, -(1/Re+gam)*minv]])
            AL_lin = LinMat(AL, [newton_a.mat.col_pardofs, 
                                 b.mat.col_pardofs]) 
            GMRes(picard.mat, res,  x=du, 
                  pre=AL_lin,
                  maxsteps=500000,
                  tol=gmr_tol, 
                  printrates=comm.rank==0)
        if pre_use == 'none':
            du.data = ns.mat.Inverse(X.FreeDofs(),
                                     inverse='pardiso') * res
        # correction for gfu
        gfu.vec.data -= dampfactor * du
    # stopping criteria 
    stopcritval = sqrt(abs(InnerProduct(du,res)))
    master_print("<A u",it,", A u",it,">_{-1}^0.5 =", stopcritval)
    if stopcritval < nonlinear_tol:
        master_print("Newton converge!")
        master_print("Newton iters",it)
        break
t1 = MPI.Wtime()
solve_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("[max_wall_time =",solve_time,"s]")

# visual the sub-domains obtained by the automatic partitioning
master_print_line("=","CHECK")
t0 = MPI.Wtime()
fesL2 = L2(mesh, order=0)
gfrank = GridFunction(fesL2)
gfrank.vec.local_vec[:] = comm.rank
div = sqrt(Integrate(div(gfu.components[0])*div(gfu.components[0]),mesh))
t1 = MPI.Wtime()
check_time = comm.allreduce(t1-t0, MPI.MAX)
master_print("divergence",div)
master_print("[max_wall_time =",check_time,"s]")

# output by pickle
output_time = 0
import os
import pickle
do_pickle = True
if do_pickle:
    output_path = os.getcwd()+"/../pickles/sns_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT PICKLE")
    t0 = MPI.Wtime()
    netgen.meshing.SetParallelPickling(True)
    pickle.dump((gfrank, gfu.components[0],gfu.components[1]),
                open(output_path+"/pickleout"+str(comm.rank), "wb"))
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",output_time,"s]")
    
# output by vtk
do_vtk = False
if do_vtk:
    output_path = os.getcwd()+"/../vtks/sns_2D"
    if rank==0 and not os.path.exists(output_path):
            os.mkdir(output_path)
    comm.Barrier() #wait until master has created the directory!!
    master_print_line("=","OUTPUT VTK")
    t0 = MPI.Wtime()
    vtk = VTKOutput(mesh, 
                    coefs=[gfrank,gfu.components[0],gfu.components[1]],
                    names=["rank","u","p"], 
                    filename=output_path+"/vtkout", subdivision=2)
    vtk.Do()
    t1 = MPI.Wtime()
    output_time = comm.allreduce(t1-t0, MPI.MAX)
    master_print("output_dirct =",output_path)
    master_print("[max_wall_time =",output_time,"s]")
    
# wall time 
wall_time = mesh_time + assemble_time + solve_time + check_time + output_time
master_print_line("=","TOTAL WALL TIME")
master_print("[tol_wall_time =",wall_time,"s]")