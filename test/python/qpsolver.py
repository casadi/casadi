from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class QPSolverTests(casadiTestCase):

  def test_OOQPsmall(self):
    try:
      OOQPSolver()
    except:
      self.message("OOQP not tested")
      return
    self.message("OOQP")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    UBA = DMatrix([2, 2, 3])
    LBA = DMatrix([-inf]*3)

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)


    solver = OOQPSolver(H.sparsity(),G.sparsity(),A.sparsity())
    solver.init()

    solver.input(QP_H).set(H)
    solver.input(QP_G).set(G)
    solver.input(QP_A).set(A)
    solver.input(QP_LBX).set(LBX)
    solver.input(QP_UBX).set(UBX)
    solver.input(QP_LBA).set(LBA)
    solver.input(QP_UBA).set(UBA)

    solver.solve()

    self.checkarray(solver.output(),DMatrix([2.0/3,4.0/3]),"OOQP")
    
if __name__ == '__main__':
    unittest.main()
