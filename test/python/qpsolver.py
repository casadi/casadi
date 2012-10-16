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
  qpsolvers.append((NLPQPSolver,{"nlp_solver": WorhpSolver, "nlp_solver_options": {}}))
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
try:
  qpsolvers.append((CplexSolver,{}))
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

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
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
      self.assertAlmostEqual(solver.output()[0],1,3,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1,3,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_COST),-6,6,str(qpsolver))
      
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,6,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([2,0,0]),str(qpsolver),digits=2)
      
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

    options = { "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
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

    options = { "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))
      if not("Cplex" in str(qpsolver)):
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

      self.assertRaises(Exception,lambda : solver.solve())

  def test_equality(self):
    self.message("Regression 452 test: equality constraints give wrong multipliers")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
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

      self.assertAlmostEqual(solver.output()[0],0.4,4,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],1.6,4,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,5,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,5,str(qpsolver))

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([3.2,0,0]),str(qpsolver),digits=5)
       
      self.assertAlmostEqual(solver.output(QP_COST)[0],-8.4,5,str(qpsolver))

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

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
      
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

      self.checkarray(solver.output(),DMatrix([5.5,5,-10]),str(qpsolver),digits=4) 
      
      self.checkarray(solver.output(QP_LAMBDA_X),DMatrix([0,0,-2.5]),str(qpsolver),digits=4)

      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([1.5]),str(qpsolver),digits=4)
       
      self.assertAlmostEqual(solver.output(QP_COST)[0],-38.375,5,str(qpsolver))
        
    
  def test_no_inequality(self):
    self.message("No inequalities present")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1]])

      


    LBA = DMatrix([0.5])
    UBA = DMatrix([0.5])

    LBX = DMatrix([-10])
    UBX = DMatrix([10])


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
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


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      if "Cplex" in str(qpsolver):
        continue
      self.message("no A: " + str(qpsolver))
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

      self.assertAlmostEqual(solver.output()[0],10,3,str(qpsolver))
      self.assertAlmostEqual(solver.output()[1],8,3,str(qpsolver))
    
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[0],0,5,str(qpsolver))
      self.assertAlmostEqual(solver.output(QP_LAMBDA_X)[1],0,5,str(qpsolver))


      self.checkarray(solver.output(QP_LAMBDA_A),DMatrix([]),str(qpsolver),digits=5)
      
      self.assertAlmostEqual(solver.output(QP_COST)[0],-34,5,str(qpsolver))
      
  def test_badscaling(self):
    return
    self.message("Badly scaled problem")
    N = 50
    H = c.diag(range(1,N+1))
    x0 = DMatrix(range(N))
    
    G = -1.0*mul(H,x0)
    print -1.0*mul(H,x0)
    print -mul(H,x0)
    A =  DMatrix(0,N)

    LBX = DMatrix([-1000]*N)
    UBX = DMatrix([1000]*N)


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
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

      solver.solve()

      self.checkarray(solver.output(),x0,str(qpsolver),digits=2)
      self.assertAlmostEqual(solver.output(QP_COST)[0],0,3,str(qpsolver))
      self.checkarray(solver.output(QP_LAMBDA_X),DMatrix.zeros(N,1),str(qpsolver),digits=4)
      
  def test_redundant(self):
    self.message("Redundant constraints")
    
    H = DMatrix([[1,-1,0],[-1,2,0],[0,0,0]])
    G = DMatrix([-2,-6,1])
    a = DMatrix([1,0,1])
    a_ = DMatrix([0,1,-2])
    
    for w0,w1 in [(0,2),(1,1),(0.1,0.6)]:
      
      A =  vertcat([a.T,a_.T,(w0*a+w1*a_).T])
        
      LBA = DMatrix([0,0,0])
      UBA = DMatrix([0.5,0.3,w0*0.5+w1*0.3])

      LBX = DMatrix([-10])
      UBX = DMatrix([10])
      
      options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
        
      for qpsolver, qp_options in qpsolvers:
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
        
        self.checkarray(solver.output(),DMatrix([-0.19230768069,1.6846153915,0.692307690769276]),str(qpsolver),digits=6)
        self.assertAlmostEqual(solver.output(QP_COST)[0],-5.850384678537,5,str(qpsolver))
        self.checkarray(solver.output(QP_LAMBDA_X),DMatrix([0,0,0]),str(qpsolver),digits=6)
        self.checkarray(mul(A.T,solver.output(QP_LAMBDA_A)),DMatrix([3.876923073076,2.4384615365384965,-1]),str(qpsolver),digits=6)
      
if __name__ == '__main__':
    unittest.main()
