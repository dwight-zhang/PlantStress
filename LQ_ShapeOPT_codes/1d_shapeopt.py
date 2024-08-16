'''1d_shapeopt.py
define:
  - float
      :SCALE: step size for shape optimization
  - list
      :MOVE_MESH: [bool: whether to move mesh, str: boundary mark for fixed boundary]
      :LIP: [bool: whether to use p-laplace equation , int: p value of p-laplace]
      :COEFF_CR: coefficient for CR equation, only relevant for LIP[0] == FALSE
      :GEO: [bool: whether to use geometric constraint,
             np.array: barycenter coordinates,
             float: volume for current mesh]
      :GEO_T: desired geometric quantaties
      :LAM: lagrange multipliers for geometric constraints
      :RHO: penalty factors for geometric constraints
      :TOL_LAM: tolerances for geometric  quantities
      :TOL_J: a sequence of tolarences for shape function J
  - function:
      :BarycenterCoord(mesh): return barycenter coordinate for mesh
author:
  ZHANG Donghang
  zdh@lsec.cc.ac.cn
'''

MOVE_MESH = [True,""]
LIP = [False,4]
SCALE = 1e-4
IT_MAX = 5000
TOL_J = np.logspace(-1,-10,10)
COEFF_CR = [20, 500]

def BarycenterCoord(mesh):
    ''' return barycenter coordinate for mesh

    :param mesh: ngsolve.comp.Mesh
    :return: numpy.array
    '''
    coord = np.zeros(mesh.dim)
    volume = Integrate(1,mesh)
    X = CoefficientFunction((x,y))
    for i in range(mesh.dim):
        coord[i] = Integrate(X[i]*dx,mesh)
        coord[i] /= volume
    return coord

GEO = [True,BarycenterCoord(MESH),Integrate(1,MESH)]
GEO_C = [True,True]
GEO_T = [np.zeros(MESH.dim),MESH_VOL]
LAM =[np.zeros(MESH.dim),0]
TOL_LAM = [1e-5, 1e-4]
RHO = [1e2, 10, 2]

def ShapeDerivative(P,Q,G,COEFF,GEO,GEO_T,LAM,RHO,mesh,BND_MARK,VEC):
    '''  DtN operator

    :param P: ngsolve.comp.GridFunction
    :param Q: ngsolve.comp.GridFunction
    :param G: ngsolve.fem.CoefficientFunction
    :param COEFF: list
    :param mesh: ngsolve.comp.Mesh
    :param BND_MARK: string
    :param VEC: ngsolve.comp.TestFunction
    :return:  ngsolve.comp.SumOfIntegrals
    '''
    # cost function
    PD = PointDefect(mesh,COEFF)
    GAM = InnerProduct(P-PD,P-PD)

    # normal derivative term
    GR_GX = G.Diff(x)
    GR_GY = G.Diff(y)
    GR_G = CoefficientFunction((GR_GX,GR_GY))

    FES_ST = MatrixValued(H1(mesh,order=2), symmetric=True)
    ST_P = GridFunction(FES_ST)
    ST_Q = GridFunction(FES_ST)
    ST_P.Interpolate(grad(P))
    ST_Q.Interpolate(grad(Q))

    N = specialcf.normal(mesh.dim)
    GAM += InnerProduct((ST_P-GR_G)*N,ST_Q*N)

    # geometric lagrange multipliers and penalty terms
    if GEO[0]:
        X = CoefficientFunction((x,y))
        if GEO_C[0]:
            DIFF_GEO = GEO[1]-GEO_T[0]
            GAM += (LAM[0][0]*X[0]+LAM[0][1]*X[1])/GEO[2]
            GAM += -np.dot(LAM[0],GEO[1])/GEO[2]
            GAM += -RHO[0]*np.dot(GEO[1],DIFF_GEO)/GEO[2]
            GAM += RHO[0]*(X[0]*DIFF_GEO[0]+X[1]*DIFF_GEO[1])/GEO[2]
        if GEO_C[1]:
            DIFF_GEO = GEO[2]-GEO_T[1]
            GAM += LAM[1]
            GAM += RHO[1]*DIFF_GEO

    return GAM*(VEC*N)*ds(definedon=mesh.Boundaries(BND_MARK))

