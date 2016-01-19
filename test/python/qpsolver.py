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
import numpy
import unittest
from types import *
from helpers import *

qpsols = []
if has_nlpsol("ipopt"):
  ipopt_options = {"fixed_variable_treatment":"relax_bounds",
                   "jac_c_constant":"yes",
                   "jac_d_constant":"yes",
                   "hessian_constant":"yes",
                   "tol":1e-12}
  qpsols.append(("nlpsol",{"nlpsol":"ipopt", "nlpsol_options.ipopt": ipopt_options},{}))

if has_qpsol("ooqp"):
  qpsols.append(("ooqp",{},{"less_digits":1}))

if has_qpsol("qpoases"):
  qpsols.append(("qpoases",{},{}))

if has_qpsol("cplex"):
  qpsols.append(("cplex",{},{}))

# if has_qpsol("sqic"):
#   qpsols.append(("sqic",{},{}))

print qpsols

class QpsolTests(casadiTestCase):

  def testboundsviol(self):
  
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf,-inf])


    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_convex: " + str(qpsol))

      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0
      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      with self.assertRaises(Exception):
        solver_out = solver(solver_in)

    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DM([-inf,5,-inf])
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_convex: " + str(qpsol))

      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      with self.assertRaises(Exception):
        solver_out = solver(solver_in)
        
  def test_scalar(self):
    # 1/2 x H x + G' x
    H = DM([1])
    G = DM([1])

    A =  DM(2)

    LBA = DM(-10)
    UBA = DM(10)

    LBX = DM([-10])
    UBX = DM([10])
    
    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_convex: " + str(qpsol))

      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],-1,max(1,6-less_digits),str(qpsol))
    
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([0]),str(qpsol),digits=max(1,6-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-0.5,max(1,6-less_digits),str(qpsol))

  def test_general_convex_dense(self):
    self.message("Convex dense QP with solvers: " + str([qpsol for qpsol,options,aux_options in qpsols]))
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_convex: " + str(qpsol))

      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(qpsol))
    
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([3+1.0/9,4.0/9,0]),str(qpsol),digits=max(1,6-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-8-2.0/9,max(1,6-less_digits),str(qpsol))
      
      solver_in["h"]=H*4

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],1,max(1,3-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,3-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["cost"],-6,max(1,6-less_digits),str(qpsol))
      
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(qpsol),digits=max(1,2-less_digits))
      
      solver_in["h"]=0
      
      if 'qcqp' in str(qpsol): continue # Singular hessian

      solver_out = solver(solver_in)
      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["cost"],-9-1.0/3,max(1,6-less_digits),str(qpsol))
      
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([10.0/3,4.0/3,0]),str(qpsol),digits=max(1,4-less_digits))

      solver_in["lba"]=[-inf]*3 #  Upper _and_ lower 
      solver_in["uba"]=[inf]*3  #  bounds infinite?

      solver_in["ubx"]=5

      if "worhp" in str(qp_options):
        with self.assertRaises(Exception):
          solver_out = solver(solver_in)
        return
      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],5,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],5,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["cost"],-40,max(1,5-less_digits),str(qpsol))
      
      self.assertAlmostEqual(solver_out["lam_x"][0],2,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],6,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([0,0,0]),str(qpsol),digits=max(1,4-less_digits))

  @memory_heavy()
  def test_general_convex_sparse(self):
    self.message("Convex sparse QP with solvers: " + str([qpsol for qpsol,options,aux_options in qpsols]))
    H = c.diag([2,1,0.2,0.7,1.3])

    H[1,2]=0.1
    H[2,1]=0.1
    
    G = DM([-2,-6,1,0,0])
    A =  DM([[1, 0,0.1,0.7,-1],[0.1, 2,-0.3,4,0.1]])
    A = sparsify(A)
    
    LBA = DM([-inf])
    UBA = DM([2, 2])

    LBX = DM([0]*5)
    UBX = DM([inf]*5)

  
    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_convex: " + str(qpsol))

      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([0.873908,0.95630465,0,0,0]),str(qpsol),digits=max(1,6-less_digits))
      
      self.checkarray(solver_out["lam_x"],DM([0,0,-0.339076,-10.0873907,-0.252185]),6,str(qpsol),digits=max(1,6-less_digits))

      self.checkarray(solver_out["lam_a"],DM([0,2.52184767]),str(qpsol),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-6.264669320767,max(1,6-less_digits),str(qpsol))

  def test_general_nonconvex_dense(self):
    self.message("Non convex dense QP with solvers: " + str([qpsol for qpsol,options,aux_options in qpsols]))
    H = DM([[1,-1],[-1,-2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])
    
    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    for qpsol, qp_options, aux_options in qpsols:
      self.message("general_nonconvex: " + str(qpsol))
      if not("cplex" in str(qpsol)):
        continue
      solver = casadi.qpsol("mysolver",qpsol, {'h':H.sparsity(),'a':A.sparsity()},qp_options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      self.assertRaises(Exception,lambda : solver(solver_in))

  def test_equality(self):
    self.message("Regression 452 test: equality constraints give wrong multipliers")
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
      
    for qpsol, qp_options, aux_options in qpsols:
      self.message("equality: " + str(qpsol))
      if "ooqp" in str(qpsol):
        continue
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':Sparsity.dense(3,2)},qp_options)

      try:
        less_digits=aux_options["less_digits"]      
      except:
        less_digits=0 

      A =  DM([[1, 1],[-1, 2],[2, 1]])

      LBA = DM([-inf]*3)
      UBA = DM([2, 2, 3])

      LBX = DM([0.5,0])
      UBX = DM([0.5,inf])

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA
      if 'worhp' in str(qp_options):
        with self.assertRaises(Exception):
          solver_out = solver(solver_in)
        return

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],1.25,max(1,6-less_digits),str(qpsol))
    
      self.assertAlmostEqual(solver_out["lam_x"][0],4.75,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([0,2,0]),str(qpsol),digits=max(1,6-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-7.4375,max(1,6-less_digits),str(qpsol))
    
      A =  DM([[1, 1],[-1, 2],[2, 1]])
      LBA = DM([2,-inf,-inf])
      UBA = DM([2, inf, inf])

      LBX = DM([-inf]*2)
      UBX = DM([inf]*2)


      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.4,max(1,4-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],1.6,max(1,4-less_digits),str(qpsol))
    
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,5-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,5-less_digits),str(qpsol))

      self.checkarray(solver_out["lam_a"],DM([3.2,0,0]),str(qpsol),digits=max(1,5-less_digits))
       
      self.assertAlmostEqual(solver_out["cost"][0],-8.4,max(1,5-less_digits),str(qpsol))

  @memory_heavy()
  def test_degenerate_hessian(self):
    self.message("Degenerate hessian")
    
    H = DM([[1,-1,0],[-1,2,0],[0,0,0]])
    H = sparsify(H)
    G = DM([-2,-6,1])
    A =  DM([[1, 1,1]])

      


    LBA = DM([0.5])
    UBA = DM([0.5])

    LBX = DM([-10])
    UBX = DM([10])

    for qpsol, qp_options, aux_options in qpsols:
      self.message("degenerate hessian: " + str(qpsol))
      if 'qcqp' in str(qpsol): continue
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0 

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([5.5,5,-10]),str(qpsol),digits=max(1,4-less_digits)) 
      
      self.checkarray(solver_out["lam_x"],DM([0,0,-2.5]),str(qpsol),digits=max(1,4-less_digits))

      self.checkarray(solver_out["lam_a"],DM([1.5]),str(qpsol),digits=max(1,4-less_digits))
       
      self.assertAlmostEqual(solver_out["cost"][0],-38.375,max(1,5-less_digits),str(qpsol))
        
    
  def test_no_inequality(self):
    self.message("No inequalities present")
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1]])

      


    LBA = DM([0.5])
    UBA = DM([0.5])

    LBX = DM([-10])
    UBX = DM([10])


    for qpsol, qp_options, aux_options in qpsols:
      self.message("no inequality: " + str(qpsol))
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)
      
      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.assertAlmostEqual(solver_out["x"][0],-0.5,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,6-less_digits),str(qpsol))
    
      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(qpsol))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(qpsol))


      self.checkarray(solver_out["lam_a"],DM([3.5]),str(qpsol),digits=max(1,6-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-3.375,max(1,6-less_digits),str(qpsol))

  def test_no_A(self):
    self.message("No A present")
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM(0,2)

    LBA = DM(0,1)
    UBA = DM(0,1)

    LBX = DM([-10])
    UBX = DM([10])


 
    for qpsol, qp_options, aux_options in qpsols:
      if "cplex" in str(qpsol):
        continue
      self.message("no A: " + str(qpsol))
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)
      
      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([10,8]),str(qpsol),digits=max(1,3-less_digits))
      
      self.checkarray(solver_out["lam_x"],DM([0,0]),str(qpsol),digits=max(1,4-less_digits))

      self.checkarray(solver_out["lam_a"],DM([]),str(qpsol),digits=max(1,5-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-34,max(1,5-less_digits),str(qpsol))
      
  def test_standard_form(self):
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([1,1]).T

    LBA = DM([-inf])
    UBA = DM([1])

    LBX = DM([-10])
    UBX = DM([10])

    for qpsol, qp_options, aux_options in qpsols:
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      
      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([-0.2,1.2]),str(qpsol),digits=max(1,3-less_digits))
      
      self.checkarray(solver_out["lam_x"],DM([0,0]),str(qpsol),digits=max(1,4-less_digits))

      self.checkarray(solver_out["lam_a"],DM([3.4]),str(qpsol),digits=max(1,5-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],-5.1,max(1,5-less_digits),str(qpsol))
  
  @memory_heavy()
  def test_badscaling(self):
    #return
    self.message("Badly scaled problem")
    N = 50
    H = c.diag(range(1,N+1))
    x0 = DM(range(N))
    
    G = -1.0*mtimes(H,x0)

    A =  DM(0,N)

    LBX = DM([-1000]*N)
    UBX = DM([1000]*N)

    for qpsol, qp_options, aux_options in qpsols:
      if 'cplex' in str(qpsol):
        continue
      if 'worhp' in str(qpsol): # works but occasionaly throws segfaults, ulimit on travis?
        continue
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],x0,str(qpsol)+str(qp_options),digits=max(1,2-less_digits))
      self.assertAlmostEqual(solver_out["cost"][0],-0.5*mtimes([x0.T,H,x0]),max(1,3-less_digits),str(qpsol))
      self.checkarray(solver_out["lam_x"],DM.zeros(N,1),str(qpsol),digits=max(1,4-less_digits))
      
  def test_redundant(self):
    self.message("Redundant constraints")
    
    H = DM([[1,-1,0],[-1,2,0],[0,0,0]])
    G = DM([-2,-6,1])
    a = DM([1,0,1])
    a_ = DM([0,1,-2])
    
    for w0,w1 in [(0,2),(1,1),(0.1,0.6)]:
      
      A =  vertcat([a.T,a_.T,(w0*a+w1*a_).T])
        
      LBA = DM([0,0,0])
      UBA = DM([0.5,0.3,w0*0.5+w1*0.3])

      LBX = DM([-10])
      UBX = DM([10])
      
      options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}
        
      for qpsol, qp_options, aux_options in qpsols:
        if 'qcqp' in str(qpsol): continue
        solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

        try:
          less_digits=aux_options["less_digits"]
        except:
          less_digits=0

        solver_in = {}        
        solver_in["h"]=H
        solver_in["g"]=G
        solver_in["a"]=A
        solver_in["lbx"]=LBX
        solver_in["ubx"]=UBX
        solver_in["lba"]=LBA
        solver_in["uba"]=UBA
        solver_out = solver(solver_in)

        self.checkarray(solver_out["x"],DM([-0.19230768069,1.6846153915,0.692307690769276]),str(qpsol),digits=max(1,6-less_digits))
        self.assertAlmostEqual(solver_out["cost"][0],-5.850384678537,max(1,5-less_digits),str(qpsol))
        self.checkarray(solver_out["lam_x"],DM([0,0,0]),str(qpsol),digits=max(1,6-less_digits))
        self.checkarray(mtimes(A.T,solver_out["lam_a"]),DM([3.876923073076,2.4384615365384965,-1]),str(qpsol),digits=max(1,6-less_digits))
        
  def test_linear(self):
    H = DM(2,2)
    A = DM([ [-1,1],[1,1],[1,-2]])
    LBA = DM([ -inf, 2, -inf ])
    UBA = DM([ 1, inf, 4 ])
    LBX = DM([ -inf, 0 ])
    UBX = DM([ inf, inf ])
    G = DM([ 2, 1 ])


    for qpsol, qp_options, aux_options in qpsols:
      if 'qcqp' in str(qpsol): continue
      solver = casadi.qpsol("mysolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([0.5,1.5]),str(qpsol),digits=max(1,5-less_digits))
      self.checkarray(solver_out["lam_x"],DM([0,0]),str(qpsol),digits=max(1,5-less_digits))

      self.checkarray(solver_out["lam_a"],DM([0.5,-1.5,0]),str(qpsol),digits=max(1,5-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],2.5,max(1,5-less_digits),str(qpsol))
      
  def test_linear2(self):
    H = DM(2,2)
    A = DM([[-1,1],[1,1],[1,-2]])
    LBA = DM([ -inf, 2, -inf ])
    UBA = DM([ 1, inf, 4 ])
    LBX = DM([ -inf, 3 ])
    UBX = DM([ inf, 3 ])
    G = DM([ 2.0, 1.0 ])


    for qpsol, qp_options, aux_options in qpsols:
      if 'qcqp' in str(qpsol): continue
      if 'nlp' in str(qpsol): continue
      solver = casadi.qpsol("msyolver",qpsol,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(solver_in)

      self.checkarray(solver_out["x"],DM([2,3]),str(qpsol),digits=max(1,5-less_digits))
      self.checkarray(solver_out["lam_x"],DM([0,-3]),str(qpsol),digits=max(1,5-less_digits))

      self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(qpsol),digits=max(1,5-less_digits))
      
      self.assertAlmostEqual(solver_out["cost"][0],7,max(1,5-less_digits),str(qpsol))
      
if __name__ == '__main__':
    unittest.main()
