#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
# 
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
# 
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# 
# 
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


    solver = OOQPSolver(H.sparsity(),A.sparsity())
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
      self.message("QPOases not tested")
      return
    self.message("QPOases")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    UBA = DMatrix([2, 2, 3])
    LBA = DMatrix([-inf]*3)

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)


    solver = QPOasesSolver(H.sparsity(),A.sparsity())
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

    solver = IpoptQPSolver(H.sparsity(),A.sparsity())
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
