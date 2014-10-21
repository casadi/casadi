#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
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
if NlpSolver.hasPlugin("ipopt"):
  qpsolvers.append(("nlp",{"nlp_solver":"ipopt", "nlp_solver_options": {"tol": 1e-12}}))

if NlpSolver.hasPlugin("ipopt"):
  qpsolvers.append(("nlp.ipopt",{"nlp_solver_options": {"tol": 1e-12}}))

if NlpSolver.hasPlugin("worhp"):
  qpsolvers.append(("nlp",{"nlp_solver": "worhp", "nlp_solver_options": {"TolOpti": 1e-12}}))

if NlpSolver.hasPlugin("worhp"):
  qpsolvers.append(("nlp.worhp",{"nlp_solver_options": {"TolOpti": 1e-12}}))

if QpSolver.hasPlugin("ooqp"):
  qpsolvers.append(("ooqp",{}))

if QpSolver.hasPlugin("qpoases"):
  qpsolvers.append(("qpoases",{}))

if QpSolver.hasPlugin("cplex"):
  qpsolvers.append(("cplex",{}))

if QpSolver.hasPlugin("sqic"):
  qpsolvers.append(("sqic",{}))

if SdpSolver.hasPlugin("dsdp"):
  qpsolvers.append(("qcqp",{"qcqp_solver":"socp","qcqp_solver_options": {"socp_solver": "sdp", "socp_solver_options": {"sdp_solver": "dsdp", "sdp_solver_options": {"gapTol":1e-10}} }}))

if SdpSolver.hasPlugin("dsdp"):
  qpsolvers.append(("qcqp.socp.sdp.dsdp",{"qcqp_solver_options": {"socp_solver_options": {"sdp_solver_options": {"gapTol":1e-10}} }}))

print qpsolvers

class QpSolverTests(casadiTestCase):

  def testboundsviol(self):
  
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf]*3)
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf,-inf])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      with self.assertRaises(Exception):
        solver.evaluate()

    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DMatrix([-inf,5,-inf])
    UBA = DMatrix([2, 2, 3])

    LBX = DMatrix([0]*2)
    UBX = DMatrix([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      with self.assertRaises(Exception):
        solver.evaluate()
        
  def test_scalar(self):
    # 1/2 x H x + G' x
    H = DMatrix([1])
    G = DMatrix([1])

    A =  DMatrix(2)

    LBA = DMatrix(-10)
    UBA = DMatrix(10)

    LBX = DMatrix([-10])
    UBX = DMatrix([10])
    
    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],-1,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-0.5,6,str(qpsolver))

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

      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],4.0/3,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([3+1.0/9,4.0/9,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-8-2.0/9,6,str(qpsolver))
      
      solver.setInput(H*4,"h")

      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput()[0],1,3,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],1,3,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("cost"),-6,6,str(qpsolver))
      
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([2,0,0]),str(qpsolver),digits=2)
      
      solver.setInput(0,"h")
      
      if 'qcqp' in str(qpsolver): continue # Singular hessian

      solver.evaluate()
      self.assertAlmostEqual(solver.getOutput()[0],2.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],4.0/3,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("cost"),-9-1.0/3,6,str(qpsolver))
      
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([10.0/3,4.0/3,0]),str(qpsolver),digits=4)

      solver.setInput([-inf]*3,"lba") #  Upper _and_ lower 
      solver.setInput([inf]*3,"uba")  #  bounds infinite?

      solver.setInput(5,"ubx")

      if "worhp" in str(qp_options):
        with self.assertRaises(Exception):
          solver.evaluate()
        return
      solver.evaluate()
      
      self.assertAlmostEqual(solver.getOutput()[0],5,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],5,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("cost"),-40,5,str(qpsolver))
      
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],2,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],6,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0,0,0]),str(qpsolver),digits=4)

  def test_general_convex_sparse(self):
    self.message("Convex sparse QP with solvers: " + str([qpsolver for qpsolver,options in qpsolvers]))
    H = c.diag([2,1,0.2,0.7,1.3])

    H[1,2]=0.1
    H[2,1]=0.1
    
    G = DMatrix([-2,-6,1,0,0])
    A =  DMatrix([[1, 0,0.1,0.7,-1],[0.1, 2,-0.3,4,0.1]])
    A = sparse(A)
    
    LBA = DMatrix([-inf])
    UBA = DMatrix([2, 2])

    LBX = DMatrix([0]*5)
    UBX = DMatrix([inf]*5)

    options = { "mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("general_convex: " + str(qpsolver))

      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()
      
      self.checkarray(solver.getOutput(),DMatrix([0.873908,0.95630465,0,0,0]),str(qpsolver),digits=6)
      
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0,-0.339076,-10.0873907,-0.252185]),6,str(qpsolver),digits=6)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0,2.52184767]),str(qpsolver),digits=6)

      self.assertAlmostEqual(solver.getOutput("cost")[0],-6.264669320767,6,str(qpsolver))

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
      self.message("general_nonconvex: " + str(qpsolver))
      if not("cplex" in str(qpsolver)):
        continue
      solver = QpSolver(qpsolver, qpStruct(h=H.sparsity(),a=A.sparsity()))
      
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      self.assertRaises(Exception,lambda : solver.evaluate())

  def test_equality(self):
    self.message("Regression 452 test: equality constraints give wrong multipliers")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      self.message("equality: " + str(qpsolver))
      if "ooqp" in str(qpsolver):
        continue
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=Sparsity.dense(3,2)))
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

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")
      if 'worhp' in str(qp_options):
        with self.assertRaises(Exception):
          solver.evaluate()
        return

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],1.25,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],4.75,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0,2,0]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-7.4375,6,str(qpsolver))
    
      A =  DMatrix([[1, 1],[-1, 2],[2, 1]])
      LBA = DMatrix([2,-inf,-inf])
      UBA = DMatrix([2, inf, inf])

      LBX = DMatrix([-inf]*2)
      UBX = DMatrix([inf]*2)


      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],0.4,4,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],1.6,4,str(qpsolver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,5,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,5,str(qpsolver))

      self.checkarray(solver.getOutput("lam_a"),DMatrix([3.2,0,0]),str(qpsolver),digits=5)
       
      self.assertAlmostEqual(solver.getOutput("cost")[0],-8.4,5,str(qpsolver))

  def test_degenerate_hessian(self):
    self.message("Degenerate hessian")
    
    H = DMatrix([[1,-1,0],[-1,2,0],[0,0,0]])
    H = sparse(H)
    G = DMatrix([-2,-6,1])
    A =  DMatrix([[1, 1,1]])

      


    LBA = DMatrix([0.5])
    UBA = DMatrix([0.5])

    LBX = DMatrix([-10])
    UBX = DMatrix([10])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
      
    for qpsolver, qp_options in qpsolvers:
      self.message("degenerate hessian: " + str(qpsolver))
      if 'qcqp' in str(qpsolver): continue
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
        
      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([5.5,5,-10]),str(qpsolver),digits=4) 
      
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0,-2.5]),str(qpsolver),digits=4)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([1.5]),str(qpsolver),digits=4)
       
      self.assertAlmostEqual(solver.getOutput("cost")[0],-38.375,5,str(qpsolver))
        
    
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
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.assertAlmostEqual(solver.getOutput()[0],-0.5,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput()[1],1,6,str(qpsolver))
    
      self.assertAlmostEqual(solver.getOutput("lam_x")[0],0,6,str(qpsolver))
      self.assertAlmostEqual(solver.getOutput("lam_x")[1],0,6,str(qpsolver))


      self.checkarray(solver.getOutput("lam_a"),DMatrix([3.5]),str(qpsolver),digits=6)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-3.375,6,str(qpsolver))

  def test_no_A(self):
    self.message("No A present")
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix.sparse(0,2)

    LBA = DMatrix.sparse(0,1)
    UBA = DMatrix.sparse(0,1)

    LBX = DMatrix([-10])
    UBX = DMatrix([10])


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      if "cplex" in str(qpsolver):
        continue
      self.message("no A: " + str(qpsolver))
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([10,8]),str(qpsolver),digits=3)
      
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(qpsolver),digits=4)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(qpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-34,5,str(qpsolver))
      
  def test_standard_form(self):
    H = DMatrix([[1,-1],[-1,2]])
    G = DMatrix([-2,-6])
    A =  DMatrix([1,1]).T

    LBA = DMatrix([-inf])
    UBA = DMatrix([1])

    LBX = DMatrix([-10])
    UBX = DMatrix([10])


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([-0.2,1.2]),str(qpsolver),digits=3)
      
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(qpsolver),digits=4)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([3.4]),str(qpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],-5.1,5,str(qpsolver))
      
  def test_badscaling(self):
    #return
    self.message("Badly scaled problem")
    N = 50
    H = c.diag(range(1,N+1))
    x0 = DMatrix(range(N))
    
    G = -1.0*mul(H,x0)

    A =  DMatrix.sparse(0,N)

    LBX = DMatrix([-1000]*N)
    UBX = DMatrix([1000]*N)


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      if 'cplex' in str(qpsolver):
        continue
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
           
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")

      solver.evaluate()

      self.checkarray(solver.getOutput(),x0,str(qpsolver)+str(qp_options),digits=2)
      self.assertAlmostEqual(solver.getOutput("cost")[0],-0.5*mul([x0.T,H,x0]),3,str(qpsolver))
      self.checkarray(solver.getOutput("lam_x"),DMatrix.zeros(N,1),str(qpsolver),digits=4)
      
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
        if 'qcqp' in str(qpsolver): continue
        solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
        for key, val in options.iteritems():
          if solver.hasOption(key):
             solver.setOption(key,val)
        solver.setOption(qp_options)
        solver.init()
        
        solver.setInput(H,"h")
        solver.setInput(G,"g")
        solver.setInput(A,"a")
        solver.setInput(LBX,"lbx")
        solver.setInput(UBX,"ubx")
        solver.setInput(LBA,"lba")
        solver.setInput(UBA,"uba")
        solver.evaluate()
        
        self.checkarray(solver.getOutput(),DMatrix([-0.19230768069,1.6846153915,0.692307690769276]),str(qpsolver),digits=6)
        self.assertAlmostEqual(solver.getOutput("cost")[0],-5.850384678537,5,str(qpsolver))
        self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0,0]),str(qpsolver),digits=6)
        self.checkarray(mul(A.T,solver.getOutput("lam_a")),DMatrix([3.876923073076,2.4384615365384965,-1]),str(qpsolver),digits=6)
        
  def test_linear(self):
    H = DMatrix.sparse(2,2)
    A = DMatrix([ [-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, inf, 4 ])
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, inf ])
    G = DMatrix([ 2, 1 ])


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      if 'qcqp' in str(qpsolver): continue
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()
      


      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([0.5,1.5]),str(qpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(qpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0.5,-1.5,0]),str(qpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],2.5,5,str(qpsolver))
      
  def test_linear2(self):
    H = DMatrix.sparse(2,2)
    A = DMatrix([[-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, inf, 4 ])
    LBX = DMatrix([ -inf, 3 ])
    UBX = DMatrix([ inf, 3 ])
    G = DMatrix([ 2.0, 1.0 ])


    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsolver, qp_options in qpsolvers:
      if 'qcqp' in str(qpsolver): continue
      if 'nlp' in str(qpsolver): continue
      solver = QpSolver(qpsolver,qpStruct(h=H.sparsity(),a=A.sparsity()))
      for key, val in options.iteritems():
        if solver.hasOption(key):
           solver.setOption(key,val)
      solver.setOption(qp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.evaluate()

      self.checkarray(solver.getOutput(),DMatrix([2,3]),str(qpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-3]),str(qpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([2,0,0]),str(qpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],7,5,str(qpsolver))
      
if __name__ == '__main__':
    unittest.main()
