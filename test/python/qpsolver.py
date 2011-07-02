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
    solver.setOption("mutol",1e-12)
    solver.setOption("artol",1e-12)
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
    
    
    self.checkarray(solver.output(QP_COST),DMatrix(-8-2.0/9),"OOQP")
    
    solver.input(QP_H).set(H*4)

    solver.evaluate()
    self.assertAlmostEqual(solver.output()[0],1,6)
    self.assertAlmostEqual(solver.output()[1],1,6)
    self.checkarray(solver.output(QP_COST),DMatrix(-6),"OOQP")
    
    solver.input(QP_UBA).set([-inf]*3)
    
    self.assertRaises(Exception,lambda : solver.evaluate())

  def test_QPOasessmall(self):
    try:
      QPOasesSolver()
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


    solver = QPOasesSolver(H.sparsity(),G.sparsity(),A.sparsity())
    solver.init()

    solver.input(QP_H).set(H)
    solver.input(QP_G).set(G)
    solver.input(QP_A).set(A)
    solver.input(QP_LBX).set(LBX)
    solver.input(QP_UBX).set(UBX)
    solver.input(QP_LBA).set(LBA)
    solver.input(QP_UBA).set(UBA)

    solver.solve()

    self.checkarray(solver.output(),DMatrix([2.0/3,4.0/3]),"qpOASES")
    
    
    self.checkarray(solver.output(QP_COST),DMatrix(-8-2.0/9),"qpOASES")
    
    solver.input(QP_H).set(H*4)

    solver.evaluate()
    self.assertAlmostEqual(solver.output()[0],1,6)
    self.assertAlmostEqual(solver.output()[1],1,6)
    self.checkarray(solver.output(QP_COST),DMatrix(-6),"qpOASES")
    
    solver.input(QP_UBA).set([-inf]*3)
    
    self.assertRaises(Exception,lambda : solver.evaluate())
  def test_IPOPTsmall(self):
    try:
      IpoptQPSolver()
    except:
      self.message("IPOPT QP not tested")
      return
    self.message("IPOPT QP")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    UBA = DMatrix([2, 2, 3])
    LBA = DMatrix([-inf]*3)

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    solver = IpoptQPSolver(H.sparsity(),G.sparsity(),A.sparsity())
    solver.setOption("convex",True)
    solver.init()

    solver.input(QP_H).set(H)
    solver.input(QP_G).set(G)
    solver.input(QP_A).set(A)
    solver.input(QP_LBX).set(LBX)
    solver.input(QP_UBX).set(UBX)
    solver.input(QP_LBA).set(LBA)
    solver.input(QP_UBA).set(UBA)

    solver.solve()

    print solver.output()
    self.assertAlmostEqual(solver.output()[0],2.0/3,6)
    self.assertAlmostEqual(solver.output()[1],4.0/3,6)
    
    
    self.assertAlmostEqual(solver.output(QP_COST)[0],-8-2.0/9,6)
    
    solver.input(QP_H).set(H*4)

    solver.evaluate()
    self.assertAlmostEqual(solver.output()[0],1,4)
    self.assertAlmostEqual(solver.output()[1],1,4)
    self.assertAlmostEqual(solver.output(QP_COST),-6,6)
    
if __name__ == '__main__':
    unittest.main()
