from ngsolve import *
from ngsolve.la import ParallelMatrix, ParallelDofs, SparseMatrixd
from ngsolve.la import CreateParallelVector
from functools import reduce

class SymmetricGS(BaseMatrix):
    def __init__ (self, smoother):
        super(SymmetricGS, self).__init__()
        self.smoother = smoother
    def Mult (self, x, y):
        y[:] = 0.0
        self.smoother.Smooth(y, x)
        self.smoother.SmoothBack(y,x)
    def Height (self):
        return self.smoother.height
    def Width (self):
        return self.smoother.height

# Okay, this is a bit annoying:
# This class is a wrapper around a block-operator,
# that simply linearizes the blockvectors and itself
# works on standard ParallelVectors.
# Hard coded for 2x2 blocks, and does currently not
# work on multidim vectors.
#
# This is needed for the (C++-side) GMRes
#
# TODO: implement something like this on C++ side ASAP
#       that way we can also get rid of unnecessary vector
#       copying
#
#    (at the very least put a ParallelFlatVector(loc_vec, pardofs)
#     into the python interface)
#
#        -- rant over ---
#
class LinMat(BaseMatrix):
    def __init__ (self, A, pd_array):
        super(LinMat, self).__init__()
        self.A = A
        comm = MPI_Init()
        dist_procs = reduce(lambda x,y:x+y, [[list(p.Dof2Proc(k)) for k in range(p.ndoflocal)] for p in pd_array])
        self.long_pardofs = ParallelDofs(dist_procs, comm)
        self.n1 = pd_array[0].ndoflocal
        self.n2 = self.n1+pd_array[1].ndoflocal

        self.x = self.A.CreateRowVector()
        self.y = self.A.CreateRowVector()

    def IsComplex(self):
        return False
    def Height(self):
        return self.n2
    def Width(self):
        return self.n2
    def CreateRowVector(self):
        return CreateParallelVector(self.long_pardofs, self.x.GetParallelStatus())
    def CreateColVector(self):
        return CreateParallelVector(self.long_pardofs, self.x.GetParallelStatus())
    def Mult(self, x, y):
        self.x[0].local_vec.data = x.local_vec[0:self.n1]
        self.x[1].local_vec.data = x.local_vec[self.n1:self.n2]
        self.x[0].SetParallelStatus(x.GetParallelStatus())
        self.x[1].SetParallelStatus(x.GetParallelStatus())
        self.y.data = self.A * self.x
        self.y[0].Cumulate()
        self.y[1].Cumulate()
        y.local_vec[0:self.n1] = self.y[0].local_vec
        y.local_vec[self.n1:self.n2] = self.y[1].local_vec
