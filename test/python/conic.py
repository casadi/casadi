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

conics = []
if has_nlpsol("ipopt"):
  ipopt_options = {"fixed_variable_treatment":"relax_bounds",
                   "jac_c_constant":"yes",
                   "jac_d_constant":"yes",
                   "hessian_constant":"yes",
                   "tol":1e-12}
  conics.append(("nlpsol",{"nlpsol":"ipopt", "nlpsol_options.ipopt": ipopt_options},{"quadratic": True}))

if has_nlpsol("worhp"):
  worhp_options = {"TolOpti":1e-13}
  conics.append(("nlpsol",{"nlpsol":"worhp", "nlpsol_options.worhp": worhp_options},{"less_digits":1,"quadratic": True}))


if has_conic("ooqp"):
  conics.append(("ooqp",{},{"less_digits":1,"quadratic": True}))

if has_conic("qpoases"):
  conics.append(("qpoases",{},{"quadratic": True}))

if has_conic("cplex"):
  conics.append(("cplex",{},{"quadratic": True}))

# if has_conic("sqic"):
#   conics.append(("sqic",{},{}))

if has_conic("clp"):
  conics.append(("clp",{"verbose":True},{"quadratic": False}))

if has_conic("activeset"):
  conics.append(("activeset",dict(max_iter=100),{"quadratic": True}))

print(conics)

class ConicTests(casadiTestCase):

  def test_general_convex_dense(self):
    self.message("Convex dense QP with solvers: " + str([conic for conic,options,aux_options in conics]))
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])

    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      self.message("general_convex: " + str(conic))

      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)
      try:
          if conic!="activeset": self.assertTrue(solver.stats()["success"])
      except:
          raise Exception(str(conic))

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([3+1.0/9,4.0/9,0]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8-2.0/9,max(1,6-less_digits),str(conic))

      solver_in["h"]=H*4

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-6,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(conic),digits=max(1,2-less_digits))

      solver_in["h"]=0

      if 'qcqp' in str(conic): continue # Singular hessian
      if 'activeset' in str(conic): continue # Singular hessian

      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-9-1.0/3,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([10.0/3,4.0/3,0]),str(conic),digits=max(1,4-less_digits))

      solver_in["lba"]=[-inf]*3 #  Upper _and_ lower
      solver_in["uba"]=[inf]*3  #  bounds infinite?

      solver_in["ubx"]=5

      if "worhp" in str(qp_options):
        with self.assertRaises(Exception):
          solver_out = solver(solver_in)
        continue
      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-40,max(1,5-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],2,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],6,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([0,0,0]),str(conic),digits=max(1,4-less_digits))

  @memory_heavy()
  def test_general_convex_sparse(self):
    self.message("Convex sparse QP with solvers: " + str([conic for conic,options,aux_options in conics]))
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


    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      self.message("general_convex: " + str(conic))

      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([0.873908,0.95630465,0,0,0]),str(conic),digits=max(1,6-less_digits))

      self.checkarray(solver_out["lam_x"],DM([0,0,-0.339076,-10.0873907,-0.252185]),6,str(conic),digits=max(1,6-less_digits))

      self.checkarray(solver_out["lam_a"],DM([0,2.52184767]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-6.264669320767,max(1,6-less_digits),str(conic))

  def test_general_nonconvex_dense(self):
    self.message("Non convex dense QP with solvers: " + str([conic for conic,options,aux_options in conics]))
    H = DM([[1,-1],[-1,-2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])

    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    for conic, qp_options, aux_options in conics:
      self.message("general_nonconvex: " + str(conic))
      if not("cplex" in str(conic)):
        continue
      solver = casadi.conic("mysolver",conic, {'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

    for conic, qp_options, aux_options in conics:
      self.message("equality: " + str(conic))
      if not aux_options["quadratic"]: continue
      if "ooqp" in str(conic):
        continue
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':Sparsity.dense(3,2)},qp_options)

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
        continue

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1.25,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],4.75,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([0,2,0]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-7.4375,max(1,6-less_digits),str(conic))

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

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.4,max(1,4-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1.6,max(1,4-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,5-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,5-less_digits),str(conic))

      self.checkarray(solver_out["lam_a"],DM([3.2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8.4,max(1,5-less_digits),str(conic))

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

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      self.message("degenerate hessian: " + str(conic))
      if 'qcqp' in str(conic): continue
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([5.5,5,-10]),str(conic),digits=max(1,4-less_digits))

      if "worhp" not in str(qp_options): self.checkarray(solver_out["lam_x"],DM([0,0,-2.5]),str(conic),digits=max(1,4-less_digits))

      if "worhp" not in str(qp_options): self.checkarray(solver_out["lam_a"],DM([1.5]),str(conic),digits=max(1,4-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-38.375,max(1,5-less_digits),str(conic))


  def test_no_inequality(self):
    self.message("No inequalities present")
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1]])




    LBA = DM([0.5])
    UBA = DM([0.5])

    LBX = DM([-10])
    UBX = DM([10])


    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      self.message("no inequality: " + str(conic))
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],-0.5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,6-less_digits),str(conic))

      self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))


      self.checkarray(solver_out["lam_a"],DM([3.5]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-3.375,max(1,6-less_digits),str(conic))

  def test_no_A(self):
    self.message("No A present")
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM(0,2)

    LBA = DM(0,1)
    UBA = DM(0,1)

    LBX = DM([-10])
    UBX = DM([10])



    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      if "cplex" in str(conic):
        continue
      self.message("no A: " + str(conic))
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([10,8]),str(conic),digits=max(1,3-less_digits))

      self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,4-less_digits))

      self.checkarray(solver_out["lam_a"],DM([]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-34,max(1,5-less_digits),str(conic))

  def test_standard_form(self):
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([1,1]).T

    LBA = DM([-inf])
    UBA = DM([1])

    LBX = DM([-10])
    UBX = DM([10])

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)


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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([-0.2,1.2]),str(conic),digits=max(1,3-less_digits))

      self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,4-less_digits))

      self.checkarray(solver_out["lam_a"],DM([3.4]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-5.1,max(1,5-less_digits),str(conic))

  @memory_heavy()
  def test_badscaling(self):
    #return
    self.message("Badly scaled problem")
    N = 50
    H = c.diag(list(range(1,N+1)))
    x0 = DM(list(range(N)))

    G = -1.0*mtimes(H,x0)

    A =  DM(0,N)

    LBX = DM([-1000]*N)
    UBX = DM([1000]*N)

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      if 'cplex' in str(conic):
        continue
      if 'worhp' in str(conic): # works but occasionaly throws segfaults, ulimit on travis?
        continue
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],x0,str(conic)+str(qp_options),digits=max(1,2-less_digits))
      self.assertAlmostEqual(solver_out["cost"][0],-0.5*mtimes([x0.T,H,x0]),max(1,3-less_digits),str(conic))
      self.checkarray(solver_out["lam_x"],DM.zeros(N,1),str(conic),digits=max(1,4-less_digits))

  def test_redundant(self):
    self.message("Redundant constraints")

    H = DM([[1,-1,0],[-1,2,0],[0,0,0]])
    G = DM([-2,-6,1])
    a = DM([1,0,1])
    a_ = DM([0,1,-2])

    for w0,w1 in [(0,2),(1,1),(0.1,0.6)]:

      A =  vertcat(*[a.T,a_.T,(w0*a+w1*a_).T])

      LBA = DM([0,0,0])
      UBA = DM([0.5,0.3,w0*0.5+w1*0.3])

      LBX = DM([-10])
      UBX = DM([10])

      options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

      for conic, qp_options, aux_options in conics:
        if not aux_options["quadratic"]: continue
        if 'qcqp' in str(conic): continue
        solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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
        solver_out = solver(**solver_in)

        self.checkarray(solver_out["x"],DM([-0.19230768069,1.6846153915,0.692307690769276]),str(conic),digits=max(1,6-less_digits))
        self.assertAlmostEqual(solver_out["cost"][0],-5.850384678537,max(1,5-less_digits),str(conic))
        self.checkarray(solver_out["lam_x"],DM([0,0,0]),str(conic),digits=max(1,6-less_digits))
        self.checkarray(mtimes(A.T,solver_out["lam_a"]),DM([3.876923073076,2.4384615365384965,-1]),str(conic),digits=max(1,6-less_digits))

  def test_linear(self):
    H = DM(2,2)
    A = DM([ [-1,1],[1,1],[1,-2]])
    LBA = DM([ -inf, 2, -inf ])
    UBA = DM([ 1, inf, 4 ])
    LBX = DM([ -inf, 0 ])
    UBX = DM([ inf, inf ])
    G = DM([ 2, 1 ])


    for conic, qp_options, aux_options in conics:
      if 'qcqp' in str(conic): continue
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([0.5,1.5]),str(conic),digits=max(1,5-less_digits))
      self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,5-less_digits))

      self.checkarray(solver_out["lam_a"],DM([0.5,-1.5,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],2.5,max(1,5-less_digits),str(conic))

  def test_linear2(self):
    H = DM(2,2)
    A = DM([[-1,1],[1,1],[1,-2]])
    LBA = DM([ -inf, 2, -inf ])
    UBA = DM([ 1, inf, 4 ])
    LBX = DM([ -inf, 3 ])
    UBX = DM([ inf, 3 ])
    G = DM([ 2.0, 1.0 ])


    for conic, qp_options, aux_options in conics:
      if 'qcqp' in str(conic): continue
      if 'nlp' in str(conic): continue
      if 'activeset' in str(conic): continue
      solver = casadi.conic("msyolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([2,3]),str(conic),digits=max(1,5-less_digits))
      self.checkarray(solver_out["lam_x"],DM([0,-3]),str(conic),digits=max(1,5-less_digits))

      self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],7,max(1,5-less_digits),str(conic))

  @requires_conic("hpmpc")
  @requires_conic("qpoases")
  def test_hpmpc(self):

    inf = 100
    T = 10. # Time horizon
    N = 4 # number of control intervals

    # Declare model variables
    x1 = MX.sym('x1')
    x2 = MX.sym('x2')
    x = vertcat(x1, x2)
    u = MX.sym('u')

    # Model equations
    xdot = vertcat(0.6*x1 - 1.11*x2 + 0.3*u-0.03, 0.7*x1+0.01)

    # Objective term
    L = x1**2 + 3*x2**2 + 7*u**2 -0.4*x1*x2-0.3*x1*u+u -x1-2*x2

    # Fixed step Runge-Kutta 4 integrator
    F = Function('F', [x, u], [x+xdot, L])

    J = F.jacobian_old(0, 0)
    # Start with an empty NLP
    w=[]
    w0 = []
    lbw = []
    ubw = []
    J = 0
    g=[]
    lbg = []
    ubg = []

    Xs = SX.sym('X', 2, 1, N+1)
    Us = SX.sym('U', 1, 1, N+1)

    for k in range(N):

        w += [Xs[k]]

        if k==0:
          lbw += [-inf, 1]
          ubw += [inf, 1]
          w0 += [0.3, 0.7]
        elif k==2:
          lbw += [0, -inf]
          ubw += [0, inf]
          w0 += [0, 0.3]
        else:
          lbw += [-inf, -inf]
          ubw += [  inf,  inf]
          w0  += [0, 1]

        w += [Us[k]]
        lbw += [-inf]
        ubw += [inf]
        w0  += [0]

        xplus, l = F(Xs[k],Us[k])
        J+= l
        # Add equality constraint
        g   += [+3*(xplus-Xs[k+1])]
        lbg += [0, 0]
        ubg += [0, 0]
        g   += [0.1*Xs[k][1]-0.05*Us[k]]
        lbg += [-0.5*k-0.1]
        ubg += [2]
    g   += [0.1*Xs[-1][1]]
    lbg += [0.1]
    ubg += [2]

    J+= mtimes(Xs[-1].T,Xs[-1])

    w += [Xs[-1]]
    lbw += [-inf, -inf]
    ubw += [  inf,  inf]
    w0  += [0, 1]

    # Create an NLP solver
    prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}


    J = Function("J",[prob["x"]],[jacobian(prob["g"],prob["x"])])
    J(w0).print_dense()


    solver_ref = qpsol('solver', 'qpoases', prob)
    solver = qpsol('solver', 'hpmpc', prob,{"tol":1e-12,"mu0":2,"max_iter":20})
    #solver = qpsol('solver', 'hpmpc', prob,{"N":N,"nx":[2]*(N+1),"nu":[1]*N,"ng":[1]*(N+1),"tol":1e-12,"mu0":2,"max_iter":20})

    sol_ref = solver_ref(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    self.checkarray(sol_ref["x"], sol["x"])
    self.checkarray(sol_ref["lam_g"], sol["lam_g"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)
    self.checkarray(sol_ref["f"], sol["f"])

    solver = nlpsol('solver', 'sqpmethod', prob,{"qpsol": "hpmpc", "qpsol_options": {"tol":1e-12,"mu0":2,"max_iter":20}})
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    self.checkarray(sol_ref["x"], sol["x"])
    self.checkarray(sol_ref["lam_g"], sol["lam_g"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)
    self.checkarray(sol_ref["f"], sol["f"])

  @requires_conic("hpmpc")
  @requires_conic("qpoases")
  def test_hpmc_timevarying(self):

    def mat(a):
      def fl(a):
        return float(a) if len(a)>0 else 0
      return sparsify(DM([list(map(fl,i.split("\t"))) for i in a.split("\n") if len(i)>0]))
    def vec(a):
      return DM(list(map(float,a.split("\n"))))
    N = 2
    A = """1	0.2	1	-1	0	0	0	0	0	0	0	0
-0.1	0.4	0	0	-1	0	0	0	0	0	0	0
0.3	0.2	0	0	0	-1	0	0	0	0	0	0
2	0	0.3	0	0	0	0	0	0	0	0	0
1	1	0.4	0	0	0	0	0	0	0	0	0
	0	0	1	4	2	1	0.3	-1	0	0	0
	0	0	3	1	0	1	0.2	0	-1	0	0
	0	0	1	1	1	1	1	0	0	0	0
	0	0	0	0	0	0	0	2	4	0	-1
	0	0	0	0	0	0	0	2	3	1	0
	0	0	0	0	0	0	0	0	0	0	3"""
    A = mat(A)
    nx = [2,3,2,1]
    nu = [1, 2,1]
    ng = [2, 1, 1, 1]
    N = 3
    H = """7	0	0.2	0	0	0	0	0	0	0	0	0
	7	0.3	0	0	0	0	0	0	0	0	0
0.2	0.3	1	0	0	0	0	0	0	0	0	0
	0	0	3	0	0	0	1	0	0	0	0
	0	0	0	2	0.1	0	0.7	0	0	0	0
	0	0	0	0.1	1	0	1	0	0	0	0
	0	0	0	0	0	1	0.1	0	0	0	0
	0	0	1	0.7	1	0.1	2	0	0	0	0
	0	0	0	0	0	0	0	6	0	1	0
	0	0	0	0	0	0	0	0	6	0	0
	0	0	0	0	0	0	0	1	0	4	0
	0	0	0	0	0	0	0	0	0	0	9
"""

    H = mat(H)
    #solver = conic('solver', 'hpmpc', {"a": A.sparsity(), "h": H.sparsity()},{"N":N,"nx":nx,"nu":nu,"ng":ng,"tol":1e-12,"mu0":2,"max_iter":20})
    solver = conic('solver', 'hpmpc', {"a": A.sparsity(), "h": H.sparsity()},{"tol":1e-12,"mu0":2,"max_iter":20})
    solver_ref = conic('solver', 'qpoases', {"a": A.sparsity(), "h": H.sparsity()})

    g = vec("""1
1
0.2
0.4
1
0.5
0.3
1
0.6
1
1
0.7""")
    lbg = vec("""0
    0
    0
    -2
    -2
    0
    0
    -2
    0
    -2
    -2""")

    ubg = vec("""0
    0
    0
    2
    2
    0
    0
    2
    0
    2
    2""")

    lbx = vec("""0.5
    0.2
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1
    -1""")
    ubx = vec("""0.5
    0.2
    1
    1
    1
    1
    1
    1
    1
    1
    1
    1""")

    sol = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx)
    sol_ref = solver_ref(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx)

    self.checkarray(sol_ref["x"], sol["x"],digits=7)
    self.checkarray(sol_ref["lam_a"], sol["lam_a"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)

    solver = conic('solver', 'hpmpc', {"a": A.sparsity(), "h": H.sparsity()},{"tol":1e-12,"mu0":2,"max_iter":20,"warm_start":True})
    sol = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"])

    self.checkarray(sol_ref["x"], sol["x"],digits=7)
    self.checkarray(sol_ref["lam_a"], sol["lam_a"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)

if __name__ == '__main__':
    unittest.main()
