'''1a_geo.py
define:
  - float:
      :ELL_HLA: half long axis for ellipse
      :ELL_HSA: half short axis for ellipse
      :MESH_VOL: volume for ellipse 
      :G_SIZE: global mesh size
  - list:
      :L_SIZE: local points' mesh size
      :POINTS: list of points for geometric description
      :BNDS: list of boundaries for domain
  - str:
      :DB_MARK: dirichlet boundary mark
      :NB_MARK: Neumman boundary mark
  - class:
       :GEO: netgen.geom2d.SplineGeometry to descripe domain
       :MESH: ngsolve.comp.Mesh to generate mesh for domain
author:
  ZHANG Donghang
  zdh@lsec.cc.ac.cn
'''
ELL_HLA = 4.0
ELL_HSA = 1.0
G_SIZE = 0.1
L_SIZE = [0.05, 0.05]
POINTS = [(ELL_HLA,0,L_SIZE[0]),
      (ELL_HLA,ELL_HSA),
      (0,ELL_HSA,L_SIZE[1]),
      (-ELL_HLA,ELL_HSA),
      (-ELL_HLA,0,L_SIZE[0]),
      (-ELL_HLA,-ELL_HSA),
      (0,-ELL_HSA,L_SIZE[1]),
      (ELL_HLA,-ELL_HSA)]
GEO = SplineGeometry()
P1,P2,P3,P4,P5,P6,P7,P8 = [GEO.AppendPoint(*P) for P in POINTS]
BNDS = [[["spline3",P1,P2,P3],"BND1"],
        [["spline3",P3,P4,P5],"BND2"],
        [["spline3",P5,P6,P7],"BND3"],
        [["spline3",P7,P8,P1],"BND4"]]
[GEO.Append(c,bc=bc,leftdomain=1,rightdomain=0) for c,bc in BNDS]
GEO.AddCircle(c=(0,0),r=0.8*ELL_HSA,bc="circle",leftdomain=2,rightdomain=1)
MESH_VOL = 2*ELL_HSA**2*pi
GEO.SetMaterial(1,"unknown")
GEO.SetMaterial(2,"known")
GEO.SetDomainMaxH(2,G_SIZE)
ELL = Mesh(GEO.GenerateMesh(maxh=G_SIZE))
MESH = ELL
DB_MARK= "BND1|BND2|BND3|BND4"
NB_MARK= ""

