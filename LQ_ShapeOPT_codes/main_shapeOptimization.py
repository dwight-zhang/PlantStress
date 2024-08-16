''' main_shapeOptimization.py
shape optimization for 2D nematic equlibria

author:
  ZHANG Donghang
  zdh@lsec.cc.ac.cn
'''

#---------------- import netgen/ngsolve ------------------
import netgen.gui
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve.internal import visoptions
from ngsolve.solvers import *

#---------------- import other module ------------------
import sys
import time
import numpy as np
from math import pi 

#---------------- deformation setting ------------------
ngsglobals.msg_level = 1
visoptions.vecfunction='GFS'
# SetVisualization(deformation=True)

#---------------- different examples ------------------
EXAMPLE = 1
# 1: elastic energy optimization
# 2: elastic and thermal energy optimization
# 3: elastic and general volume energy optimization
if EXAMPLE == 1:
    exec(open("1a_geo.py").read())
    exec(open("1b_pde.py").read())
    exec(open("1c_cost.py").read())
    exec(open("1d_shapeopt.py").read())
elif EXAMPLE == 2:
    exec(open("2a_geo.py").read())
    exec(open("2b_pde.py").read())
    exec(open("2c_cost.py").read())
    exec(open("2d_shapeopt.py").read())
elif EXAMPLE == 3:
    exec(open("3a_geo.py").read())
    exec(open("3b_pde.py").read())
    exec(open("3c_cost.py").read())
    exec(open("3d_shapeopt.py").read())
elif EXAMPLE == 4:
    print("Please use other main file.")
    sys.exit(0)
else:
    print("Invalid example number. Insert integer between 1 and 7")
    sys.exit()

#---------------- move mesh function ------------------
def MoveNGmesh(DISP, mesh):
    ''' move mesh by displace vector

    :parm DISP: ngsolve.comp.GridFunction
    :param mesh: ngsolve.comp.Mesh
    '''
    for p in mesh.ngmesh.Points():
        mip = mesh(p[0],p[1])
        v = DISP(mip)
        p[0] += v[0]
        p[1] += v[1]
    mesh.ngmesh.Update()

#---------------- checking input ------------------
def checkInput(stream=None):
    if stream is None:
        stream = sys.stdout
    #mesh information
    stream.write("Nverts == "+str(MESH.nv)+"\n")
    stream.write("Nelems == "+str(MESH.ne)+"\n")
    stream.write("Nfaces == "+str(MESH.nface) +"\n")
    stream.write("DB_MARK == "+DB_MARK+"\n")
    stream.write("NB_MARK == "+NB_MARK+"\n")
    #FES information
    stream.write("COEFF == "+str(COEFF_ENERGY[0])+" , "+str(COEFF_ENERGY[1])+"\n")
    stream.write("ORDER == "+str(ORDER)+"\n")
    stream.write("Ndofs == "+str(FES.ndof)+"\n")
    #optimization parameters
    stream.write("SCALE  == "+str(SCALE)+"\n")
    stream.write("LAPLACE == "+str(LIP[0])+"\n")
    if LIP[0]:
        stream.write("    p == "+str(LIP[1])+"\n")
    else:
        stream.write("    p == "+str(2)+"\n")
        stream.write("    alpha_CR == "+str(COEFF_CR[0])+"\n")
        stream.write("    gamma_CR == "+str(COEFF_CR[1])+"\n")
    stream.write("MOVE_MESH == "+str(MOVE_MESH[0])+"\n")
    stream.write("MOVE_BND == "+str(DB_MARK)+" - "+str(MOVE_MESH[1])+"\n")
    stream.write("GEOMETRIC == "+str(GEO[0])+"\n")
    if GEO[0]:
        if GEO_C[0]:
            stream.write("    Constraints ==  Barycenter\n")
            stream.write("    Desired barycenter  == "+str(GEO_T[0])+"\n")
            stream.write("    Penalty == "+str(RHO[0])+"\n")
            stream.write("    Increase == "+str(RHO[2])+"\n")
            stream.write("    tolerance == "+str(TOL_LAM[0])+"\n")
        if GEO_C[1]:
            stream.write("    Constraints ==  Volume \n")
            stream.write("    Desired volume == "+str(GEO_T[1])+"\n")
            stream.write("    Penalty == "+str(RHO[1])+"\n")
            stream.write("    Increase == "+str(RHO[2])+"\n")
            stream.write("    tolerance == "+str(TOL_LAM[1])+"\n")

#---------------- write parameters ------------------
FOLDER = "output"
FILE = FOLDER + "/"+str(EXAMPLE)+"_PARAM.txt"
OWN = open(FILE,"w+")
checkInput(OWN)
OWN.close()
#writes parameters to console.
checkInput()

#---------------- vector field ------------------
VEC = VectorH1(MESH, order=2, dirichlet = MOVE_MESH[1])
FREE = BitArray(VEC.ndof+2*FES.ndof)
for i in range(VEC.ndof):
    FREE[i] = not VEC.FreeDofs()[i]
for i in range(FES.ndof):
    FREE[VEC.ndof+i] = FES.FreeDofs()[i]
    FREE[VEC.ndof+FES.ndof+i] = FES.FreeDofs()[i]

#---------------- Hilbert Extension Regularition (HER) ------------------
PHIVEC = VEC.TrialFunction()
PSIVEC = VEC.TestFunction()
HER = BilinearForm(VEC)
if LIP[0] == True:
    HER += InnerProduct(grad(PHIVEC),grad(PSIVEC))*InnerProduct(grad(PHIVEC),grad(PSIVEC))*dx
else:
    HER += COEFF_CR[0]*InnerProduct(grad(PHIVEC),grad(PSIVEC))*dx+InnerProduct(PHIVEC,PSIVEC)*dx
    HER += COEFF_CR[1]*(grad(PHIVEC)[0,0]-grad(PHIVEC)[1,1])*(grad(PSIVEC)[0,0]-grad(PSIVEC)[1,1])*dx
    HER += COEFF_CR[1]*(grad(PHIVEC)[1,0]+grad(PHIVEC)[0,1])*(grad(PSIVEC)[1,0]+grad(PSIVEC)[0,1])*dx

#---------------- grid functions for deformation space ------------------
GFS = GridFunction(VEC)
GFX = GridFunction(VEC)
if MESH.dim==2:
    GFS.Set((0,0))
elif mesh.dim == 3:
    GFS.Set((0,0,0))

#---------------- grid functions for state space ------------------
U,V = FES.TnT()
GFP = GridFunction(FES)
GFQ = GridFunction(FES)
rLdG= BilinearForm(FES)
rLdG += EnergyVolumeIntegral(U,V,COEFF_ENERGY)
rLdG += EnergyBoundaryIntegral(U,V,COEFF_ENERGY)
PD = PointDefect(MESH,COEFF_ENERGY)
DIFF_LSD = LinearForm(FES)
DIFF_LSD += CostVolumeIntegral(GFP,PD).Diff(GFP,FES.TestFunction())
DIFF_LSD += CostBoundaryIntegral(GFP,PD).Diff(GFP,FES.TestFunction())

#---------------- prepare output ------------------
JV_FILE = FOLDER+"/"+str(EXAMPLE)+"_JVAL.m"
OWN_JV = open(JV_FILE,"w+")
OWN_JV.write("Jvals = [")
OWN_JV.close()

input("press enter to start optimization")
TIME = time.time()
with TaskManager():
    # initial 
    GFP_DB_FUNC = UniaxialBoundary(MESH,COEFF_ENERGY)
    GFQ_DB_FUNC = CoefficientFunction((0,0))
    GFP = SolveStatePDE(FES,rLdG,DB_MARK,GFP_DB_FUNC)
    GFQ = SolveAdjointPDE(FES,rLdG,DIFF_LSD,DB_MARK,GFQ_DB_FUNC)
    # VTK = VTKOutput(MESH,coefs=[GFP,GFQ],names=["STATE","ADJOINT"],
                # filename="initial",subdivision=2)
    # VTK.Do()
    Draw(GFP,MESH,"STATE")
    Draw(GFQ,MESH,"ADJOINT")
    input("output inital solutions, press entor to cont'd")
    Jinit_ag = Integrate(CostVolumeIntegral(GFP,PD)+
                         CostBoundaryIntegral(GFP,PD), MESH)
    if GEO[0]:
        GEO[1] = BarycenterCoord(MESH)
        GEO[2] = Integrate(1,MESH)
        DIFF_GEO = [GEO[1]-GEO_T[0],GEO[2]-GEO_T[1]]
        NORM = [np.linalg.norm(DIFF_GEO[0]),abs(DIFF_GEO[1])]
        if GEO_C[0]:
            Jinit_ag += RHO[0]/2*NORM[0]**2
        if GEO_C[1]:
            Jinit_ag += RHO[1]/2*NORM[1]**2
    print("===init", ' cost', Jinit_ag.real )
    Jold_ag = Jinit_ag
    Jinit_so = Jinit_ag
    Jold_so = Jinit_so
    count = 0

    for k in range(len(TOL_J)):
        # shape optimization problem
        # object tolerance
        print("===so-tol ", '%-.6e' % TOL_J[k])
        for i in range(IT_MAX):
            # initial mesh deformation
            if MOVE_MESH[0] == False:
                MESH.SetDeformation(GFS)
            # solve state and adjoint equations
            GFP = SolveStatePDE(FES,rLdG,DB_MARK,GFP_DB_FUNC)
            GFQ = SolveAdjointPDE(FES,rLdG,DIFF_LSD,DB_MARK,GFQ_DB_FUNC)
            # get deformation vector 
            if GEO[0]:
                GEO[1] = BarycenterCoord(MESH)
                GEO[2] = Integrate(1,MESH)
            FX =LinearForm(ShapeDerivative(GFP,GFQ,GFP_DB_FUNC,COEFF_ENERGY,GEO,GEO_T,
                                       LAM,RHO,MESH,DB_MARK,PSIVEC)).Assemble()
            HER.Assemble()
            INV_HER = HER.mat.Inverse(freedofs=VEC.FreeDofs(),inverse="umfpack")
            GFX.vec.data = INV_HER * FX.vec
            # wirte file
            NLTWO_GFX = Norm(GFX.vec)
            OWN_JV = open(JV_FILE,"a+")
            OWN_JV.write(str(i+count) + " " + str(Jold_so) + " " +
                         str(NLTWO_GFX)+ "; \n" ) 
            OWN_JV.close()
            if MOVE_MESH[0] == False:
                MESH.UnsetDeformation()
            # calculate new geometric parameters and objective on deformated mesh
            GFS.vec.data = GFS.vec-SCALE*GFX.vec
            if MOVE_MESH[0] == False:
                MESH.SetDeformation(GFS)
            else:
                MoveNGmesh(-SCALE*GFX,MESH)
            Jnew_so = Integrate(CostVolumeIntegral(GFP,PD)+
                                CostBoundaryIntegral(GFP,PD),MESH)
            if GEO[0]:
                GEO[1] = BarycenterCoord(MESH)
                GEO[2] = Integrate(1,MESH)
                DIFF_GEO = [GEO[1]-GEO_T[0],GEO[2]-GEO_T[1]]
                NORM = [np.linalg.norm(DIFF_GEO[0]),abs(DIFF_GEO[1])]
                if GEO_C[0]:
                    Jnew_so += RHO[0]/2*NORM[0]**2
                if GEO_C[1]:
                    Jnew_so += RHO[1]/2*NORM[1]**2
            if MOVE_MESH[0] == False:
                MESH.UnsetDeformation()
            # Stop criteria
            R_TOL = abs(Jnew_so.real-Jold_so.real)/abs(Jinit_so.real)
            if GEO[0] == False:
                TOL = TOL_J[5]
            else:
                TOL = TOL_J[k]
            if R_TOL < TOL:
                count = i
                Draw(GFS,MESH,'DEFOEM')
                Draw(GFP,MESH,'STATE')
                Draw(GFQ,MESH,'ADJOINT')
                print("-----so-steps ", '%d'% count, ' R_TOL','%-.6e' % R_TOL, 
                      'converged with J = %-.6e' % Jnew_so.real)
                input('-----cont,d')
                break
            else:
                Jold_so = Jnew_so
                print("-----so-it ", '%d'% i, ' R_TOL','%-.6e' % R_TOL, 
                      'J = %-.6e' % Jnew_so.real, '| nabla J | = %-.6e' %
                      NLTWO_GFX) 
        # Augmented Lagrange optimization
        if GEO[0] == False:
            VTK = VTKOutput(MESH,coefs=[GFP,GFQ],names=["STATE","ADJOINT"],
                            filename="R_TOL"+str(TOL),subdivision=2)
            VTK.Do()
            print("===output solutions, break with no geometeric constraints")
            break
        if MOVE_MESH[0] == False:
            MESH.SetDeformation(GFS)
        GEO[2] = Integrate(1,MESH)
        GEO[1] = BarycenterCoord(MESH)
        DIFF_GEO = [GEO[1]-GEO_T[0],GEO[2]-GEO_T[1]]
        NORM = [np.linalg.norm(DIFF_GEO[0]),abs(DIFF_GEO[1])]
        if NORM[0] > TOL_LAM[0]:
            RHO[0] *= RHO[2]
        else:
            LAM[0] = np.add(LAM[0],RHO[0]*DIFF_GEO[0])
        if NORM[1] > TOL_LAM[1]:
            RHO[1] *= RHO[2]
        else:
            LAM[1] += RHO[1]*DIFF_GEO[1]
        # Stop criteria 
        Jnew_ag = Integrate(CostVolumeIntegral(GFP,PD)+
                             CostBoundaryIntegral(GFP,PD), MESH)
        if GEO_C[0]:
            Jnew_ag += RHO[0]/2*NORM[0]**2
        if GEO_C[1]:
            Jnew_ag += RHO[1]/2*NORM[1]**2
        if MOVE_MESH[0] == False:
            MESH.UnsetDeformation()
        # Stop criteria
        R_TOL = abs(Jnew_ag.real-Jold_ag.real)/abs(Jinit_ag.real)
        if R_TOL < TOL_J[-1]:
            print("===al-steps ", k, "converged with LAM = ", '%-.6e' %
                  NORM[0], ' , %-.6e' % NORM[1])
            Draw(GFS,MESH,'DEFORM')
            Draw(GFP,MESH,'STATE')
            Draw(GFQ,MESH,'ADJOINT')
            break
        else:
            Jold_ag = Jnew_ag
            Jinit_so = Jold_ag
            print("===al-it ", '%d' % k, ' R_TOL','%-.6e'% R_TOL,
                  'J = %-.6e' % Jnew_so.real, 'LAM = %-.6e' % NORM[0], ' , %-.6e' %
                  NORM[1])
   # wirte file 
    OWN_JV = open(JV_FILE,"a+")
    OWN_JV.write("]; \n")
    OWN_JV.close()
TIME = time.time() - TIME

#write fiel
OWN_JV = open(JV_FILE,"a+")
OWN_JV.write("yyaxis left; plot(Jvals(:,1), Jvals(:,2)); ylabel('J') \n")
OWN_JV.write("yyaxis right; plot(Jvals(:,1), Jvals(:,3),'--'); ylabel('$\|\ nabla J\|$', 'Interpreter', 'latex') \n")
OWN_JV.write("title('History');  \n xlabel('Iterations') \n")
OWN_JV.write("legend('J','$\|\ nabla J\|$','Interpreter','latex') \n")
OWN_JV.write("%  optimization time: " + str(TIME) + " seconds")
OWN_JV.close()

breakpoint = input("Press any key to quit")
