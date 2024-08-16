'''1b_pde.py
define:
  - int:
      :ORDER: degree of polynomial for picewise continuous finite element
  - list:
      :COEFF_ENERGY: coefficent in reduced Landau-de Gennes energy
  - class:
      :FES_P: ngsolve.comp.FESpace for state and adjoint variables
  - function:
      :EnergyVolumeIntegral(U,V,COEFF): volume integral of Landau-de Gennes energy
      :EnergyBoundaryIntegral(U,V,COEFF): boundary integral of Landau-de Gennes energy
      :UniaxialBoundary(N,COEFF): uniaxial tensor on boundary
      :SolveStatePDE(FES,ENERGY,BND_MARK,BND_FUNC): solve PDEs for state variable
      :SolveAdjointPDE(FES,ENERGY,DIFF_COST,BND_MARK,BND_FUNC): solve PDEs for adjoint variable
author:
  ZHANG Donghang
  zdh@lsec.cc.ac.cn
'''

ORDER = 1
COEFF_ENERGY = [5/sqrt(2), 64/35];
FES = H1(MESH, order=ORDER, dirichlet=DB_MARK, dim=MESH.dim)

def EnergyVolumeIntegral(U,V,COEFF):
    '''return volume integral in energy

    :param U: ngsolve.comp.FESpace.TrialFunction
    :param V: ngsolve.comp.FESpace.TrialFunction
    :COEFF: list
    :return:  ngsolve.comp.SumOfIntegrals
    '''
    return (InnerProduct(grad(U),grad(V))+COEFF[0]**2*0*x)*dx

def EnergyBoundaryIntegral(U,V,COEFF):
    '''return boundary integral in energy

    :param U: ngsolve.comp.FESpace.TrialFunction
    :param V: ngsolve.comp.FESpace.TrialFunction
    :COEFF: list
    :return:  ngsolve.comp.SumOfIntegrals
    '''
    return 0*COEFF[1]*x*ds

def UniaxialBoundary(mesh,COEFF):
    ''' return dirichlet uniaxial boundary condition

    remark: N.dim = ngsolve.comp.Mesh.dim

    :param n: ngsolve.fem.CoefficientFunction
    :COEFF: list
    :return: ngsolve.fem.CoefficientFunction
    '''
    N = specialcf.normal(mesh.dim)
    return CoefficientFunction((COEFF[1]*(N[0]*N[0]-0.5),COEFF[1]*N[0]*N[1]))

def SolveStatePDE(FES,ENERGY,BND_MARK,BND_FUNC):
    ''' solve PDEs, and return state variable

    :param FES: ngsolve.comp.FESpace
    :param ENERGY: ngsolve.comp.BilinearForm
    :param BND_MARK: string
    "param BND_FUNC: ngsolve.fem.CoefficientFunction
    :return: ngsolve.comp.GridFunction
    '''
    # initial gridfunction and residual vector
    P = GridFunction(FES)
    R = P.vec.CreateVector()

    # solve state equation
    ENERGY.Assemble()
    INV_MAT = ENERGY.mat.Inverse(freedofs=FES.FreeDofs(), inverse="umfpack")
    P.Set(BND_FUNC, definedon=FES.mesh.Boundaries(BND_MARK))
    R.data = -ENERGY.mat*P.vec
    P.vec.data += INV_MAT*R

    return P

def SolveAdjointPDE(FES,ENERGY,DIFF_COST,BND_MARK,BND_FUNC):
    ''' solve PDEs, and return adjoint variable

    :param FES: ngsolve.comp.FESpace
    :param ENERGY: ngsolve.comp.BilinearForm
    :param DIFF_COST: ngsolve.comp.LinearForm
    :param BND_MARK: string
    "param BND_FUNC: ngsolve.fem.CoefficientFunction
    :return: ,ngsolve.comp.GridFunction
    '''
    # initial gridfunction and residual vector
    Q = GridFunction(FES)
    R = Q.vec.CreateVector()

    #solve adjoint equation
    ENERGY.Assemble()
    INV_MAT = ENERGY.mat.Inverse(freedofs=FES.FreeDofs(), inverse="umfpack")
    DIFF_COST.Assemble();
    Q.Set(BND_FUNC, definedon=FES.mesh.Boundaries(BND_MARK))
    R.data = DIFF_COST.vec-ENERGY.mat*Q.vec
    Q.vec.data += INV_MAT*R

    return Q

