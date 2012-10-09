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

qpsolvers = []
try:
  qpsolvers.append((NLPQPSolver,{"nlp_solver":IpoptSolver, "nlp_solver_options": {"tol": 1e-12}}))
except:
  pass
try:
  #qpsolvers.append((NLPQPSolver,{"nlp_solver": WorhpSolver, "nlp_solver_options": {"TolOpti": 1e-12}}))
  pass
except:
  pass
try:
  qpsolvers.append((OOQPSolver,{}))
except:
  pass
try:
  qpsolvers.append((QPOasesSolver,{}))
except:
  pass

class QPSolverTests(casadiTestCase):

  def test_general_convex_dense(self):
    self.message("Convex dense QP with solvers: " + str([qpsolver for qpsolver,options in qpsolvers]))
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([3+1.0/9,4.0/9,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-8-2.0/9,6,str(qpsolver))
      
      solver.input(QP_H).set(H*4)

      solver.evaluate()
      self.assertAlmostEqual(solver.output()[0],1,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_COST),-6,6,str(qpsolver))
      
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([2,0,0]),str(qpsolver),digits=4)
      
      solver.input(QP_H).set(0)

      solver.evaluate()
      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_COST),-9-1.0/3,6,str(qpsolver))
      
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([10.0/3,4.0/3,0]),str(qpsolver),digits=4)

      solver.input(QP_LBA).set([-inf]*3) #  Upper _and_ lower 
      solver.input(QP_UBA).set([inf]*3)  #  bounds infinite?

      solver.input(QP_UBX).setAll(5)

      solver.evaluate()
      self.assertAlmostEqual(solver.output()[0],5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_COST),-40,5,str(qpsolver))
      
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],2,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],6,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([0,0,0]),str(qpsolver),digits=4)

  def test_general_convex_sparse(self):
    self.message("Convex sparse QP with solvers: " + str([qpsolver for qpsolver,options in qpsolvers]))
    H = c.diag([2,1,0.2,0.7,1.3])

    H[1,2]=0.1
    H[2,1]=0.1
    
    G = DMatrix([-2,-6,1,0,0])
    A =  DMatrix([[1, 0,0.1,0.7,-1],[0.1, 2,-0.3,4,0.1]])
    makeSparse(A)
    
    LBA = DMatrix([-inf])
    UBA = DMatrix([2, 2])

    LBX = DMatrix([0]*5)
    UBX = DMatrix([inf]*5)

    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()
      
      self.checkarray(solver.output(),DMatrix([0.873908,0.95630465,0,0,0]),str(qpsolver),digits=6)
      
      self.checkarray(solver.output(QP_LAMBDA_X),DMatrix([0,0,-0.339076,-10.0873907,-0.252185]),6,str(qpsolver),digits=6)

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([0,2.52184767]),str(qpsolver),digits=6)

      self.assertAlmostEqual(solver.output(QP_COST)[0],-6.264669320767,6,str(qpsolver))

  def test_general_nonconvex_dense(self):
    self.message("Non convex dense QP with solvers: " + str([qpsolver for qpsolver,options in qpsolvers]))
    H = DMatrix([[1,-1],[-1,-2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],4.0/3,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([4+8.0/9,20.0/9,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-10-16.0/9,6,str(qpsolver))

  def test_equality(self):
    self.message("Regression 452 test: equality constraints give wrong multipliers")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])

    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("equality: " + str(qpsolver))
      solver = qpsolver(H.sparsity(),sp_dense(3,2))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      
      
      A =  DMatrix([[1, 1],[-1, 2],[2, 1]])

      LBA = DMatrix([-inf]*3)
      UBA = DMatrix([2, 2, 3])

      LBX = DMatrix([0.5,0])
      UBX = DMatrix([0.5,inf])

      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.25,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],4.75,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([0,2,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-7.4375,6,str(qpsolver))
    
      A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
      LBA = DMatrix([2,-inf,-inf])
      UBA = DMatrix([2, inf, inf])

      LBX = DMatrix([-inf]*2)
      UBX = DMatrix([inf]*2)


      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],0.4,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.6,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([3.2,0,0]),str(qpsolver),digits=6)
       
      self.assertAlmostEqual(solver.output(QP_COST)[0],-8.4,6,str(qpsolver))

  def test_degenerate_hessian(self):
    self.message("Degenerate hessian")
    
    H = DMatrix([[1,-1,0],[-1,2,0],[0,0,0]])
    makeSparse(H)
    G = DMatrix([-2,-6,1])
    A =  DMatrix([[1, 1,1]])

      


    LBA = DMatrix([0.5])
    UBA = DMatrix([0.5])

    LBX = DMatrix([-10])
    UBX = DMatrix([10])

    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
      
    for qpsolver, qp_options in qpsolvers:
      self.message("degenerate hessian: " + str(qpsolver))
      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
        
      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.checkarray(solver.output(),DMatrix([5.5,5,-10]),str(qpsolver),digits=6) 
      
      self.checkarray(solver.output(QP_LAMBDA_X),DMatrix([0,0,-2.5]),str(qpsolver),digits=6)

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([1.5]),str(qpsolver),digits=6)
       
      self.assertAlmostEqual(solver.output(QP_COST)[0],-38.375,6,str(qpsolver))
        
    
  def test_no_inequality(self):
    self.message("No inequalities present")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1]])

      


    LBA = DMatrix([0.5])
    UBA = DMatrix([0.5])

    LBX = DMatrix([-10])
    UBX = DMatrix([10])


    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("no inequality: " + str(qpsolver))
      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],-0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))


      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([3.5]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-3.375,6,str(qpsolver))

  def test_no_A(self):
    self.message("No A present")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix(0,2)

    LBA = DMatrix(0,1)
    UBA = DMatrix(0,1)

    LBX = DMatrix([-10])
    UBX = DMatrix([10])


    options = {"convex": True, "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("no A: " + str(qpsolver))
      if 'NLP' in str(qpsolver):
        continue
      solver = qpsolver(H.sparsity(),A.sparsity())
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.input(QP_H).set(H)
      solver.input(QP_G).set(G)
      solver.input(QP_A).set(A)
      solver.input(QP_LBX).set(LBX)
      solver.input(QP_UBX).set(UBX)
      solver.input(QP_LBA).set(LBA)
      solver.input(QP_UBA).set(UBA)

      solver.solve()

      self.assertAlmostEqual(solver.output()[0],10,6,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],8,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))


      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-34,6,str(qpsolver))
      
if __name__ == '__main__':
    unittest.main()
