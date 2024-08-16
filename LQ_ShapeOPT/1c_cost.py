'''1c_cost.py
define:
  - function:
      :PointDefectTensor(COEFF): point defect tensor
      :CostVolumeIntegral(P,PD): volume integral in cost functional
      :CostBoundaryIntegral(P,PD): boundary integral in cost functional
author:
  ZHANG Donghang
  zdh@lsec.cc.ac.cn
'''
def PointDefect(mesh,COEFF):
    ''' return point defector symmetric tensor

    :COEFF: list
    :return: ngsolve.fem.CoefficientFunction
    '''
    cs = cos(atan(y/x))
    si = sin(atan(y/x))
    domain_values = {'known': (COEFF[1]*(cs*cs-0.5),COEFF[1]*cs*si),
                     'unknown': (0,0)}
    return mesh.MaterialCF(domain_values)

def CostVolumeIntegral(P,PD):
    ''' return voulme integral in cost functional

    :param P: ngsolve.comp.GridFunction
    :return:  ngsolve.comp.SumOfIntegrals
    '''
    return InnerProduct(P-PD,P-PD)*dx

def CostBoundaryIntegral(P,PD):
    ''' return boundary integral in cost functional

    :param P: ngsolve.comp.GridFunction
    :return:  ngsolve.comp.SumOfIntegrals
    '''
    return 0*x*ds 
    # n = specialcf.normal(mesh.dim)
    # u = gfu.components[0]
    # p = gfu.components[1]
    # fesstress = MatrixValued(H1(mesh,order=2), symmetric=True)
    # gfstress = GridFunction(fesstress)
    # gfstress.Interpolate(grad(u))
    # return InnerProduct(mu*gfstress*n - p*Id(mesh.dim)*n,d)*ds(definedon=mesh.Boundaries("cyl"))
