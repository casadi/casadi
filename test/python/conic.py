#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
import os
import platform

conics = []

if "SKIP_IPOPT_TESTS" not in os.environ and has_nlpsol("ipopt"):
  ipopt_options = {"fixed_variable_treatment":"relax_bounds",
                   "jac_c_constant":"yes",
                   "jac_d_constant":"yes",
                   "hessian_constant":"yes",
                   "tol":1e-12,
                   "print_level":0}
  conics.append(("nlpsol",{"nlpsol":"ipopt", "nlpsol_options.ipopt": ipopt_options},{"quadratic": True, "dual": True, "soc": False, "codegen": False, "discrete": False, "sos":False}))

if "SKIP_WORHP_TESTS" not in os.environ and has_nlpsol("worhp"):
  worhp_options = {"TolOpti":1e-13}
  conics.append(("nlpsol",{"nlpsol":"worhp", "nlpsol_options.worhp": worhp_options},{"less_digits":1,"quadratic": True, "dual": False, "soc": False, "codegen": False, "discrete": False, "sos":False}))

if "SKIP_OOQP_TESTS" not in os.environ and has_conic("ooqp"):
  conics.append(("ooqp",{},{"less_digits":1,"quadratic": True, "dual": True, "soc": False, "codegen": False, "discrete": False, "sos":False}))


if "SKIP_QPOASES_TESTS" not in os.environ and has_conic("qpoases"):
  conics.append(("qpoases",{"printLevel":"low"},{"quadratic": True, "dual": True, "soc": False, "codegen": False,"discrete": False, "sos":False}))


if "SKIP_CPLEX_TESTS" not in os.environ and has_conic("cplex"):
  conics.append(("cplex",{"cplex": {"CPX_PARAM_BARQCPEPCOMP": 1e-10,"CPX_PARAM_BAREPCOMP":1e-10}},{"quadratic": True, "dual": True, "soc": True, "codegen": False, "discrete": True, "sos": True}))

if "SKIP_OSQP_TESTS" not in os.environ and has_conic("osqp"):
  options = ["-Wno-unused-variable"]
  if os.name=='nt':
    options = []
  codegen = {"extralibs": ["osqp"], "extra_options": options, "std": "c99"}
  conics.append(("osqp",{"osqp":{"alpha":1,"eps_abs":1e-8,"eps_rel":1e-8}},{"quadratic": True, "dual": True, "codegen": codegen,"soc":False,"discrete":False}))

if "SKIP_SUPERSCS_TESTS" not in os.environ and has_conic("superscs"):
  conics.append(("superscs",{"superscs": {"eps":1e-9,"do_super_scs":1, "verbose":0}},{"quadratic": True, "dual": False, "codegen": False,"soc":True,"discrete":False}))

# No solution for licensing on travis

if "SKIP_GUROBI_TESTS" not in os.environ and has_conic("gurobi"):
  conics.append(("gurobi",{"gurobi": {"BarQCPConvTol":1e-9}},{"quadratic": True, "less_digits": True, "dual": True, "soc": True, "codegen": False,"discrete":True,"sos":True}))

# if has_conic("sqic"):
#   conics.append(("sqic",{},{}))

if "SKIP_CLP_TESTS" not in os.environ and has_conic("clp"):
  conics.append(("clp",{"verbose":True},{"quadratic": False, "dual": True, "soc": False, "codegen": False, "discrete": False, "sos":False}))

if "SKIP_CBC_TESTS" not in os.environ and has_conic("cbc"):
  conics.append(("cbc",{"verbose":True},{"quadratic": False, "dual": True, "soc": False, "codegen": False, "discrete": True, "sos":True}))

if "SKIP_QRQP_TESTS" not in os.environ and has_conic("qrqp"):
  codegen = {"std":"c99"}
  conics.append(("qrqp",{"max_iter":20,"print_header":False,"print_iter":False},{"quadratic": True, "dual": True, "soc": False, "codegen": codegen, "discrete": False, "sos":False}))

if "SKIP_PROXQP_TESTS" not in os.environ and has_conic("proxqp"):
  conics.append(("proxqp",{"proxqp":{"eps_abs":1e-11,"max_iter":1e4, "backend": "sparse"}}, {"quadratic": True, "dual": True, "soc": False, "codegen": False,"discrete":False,"sos":False}))
  conics.append(("proxqp",{"proxqp":{"eps_abs":1e-11,"max_iter":1e4, "backend": "dense"}}, {"quadratic": True, "dual": True, "soc": False, "codegen": False,"discrete":False,"sos":False}))

if "SKIP_QPALM_TESTS" not in os.environ and has_conic("qpalm"):
  eps = 1e-8
  conics.append(("qpalm",{"qpalm":{"eps_abs":eps,"eps_rel":eps,"eps_abs_in":eps,"eps_prim_inf":eps}},{"quadratic": True, "dual": True, "soc": False, "codegen": False, "discrete": False, "sos":False}))

if "SKIP_HIGHS_TESTS" not in os.environ and has_conic("highs"):
    codegen = {"extralibs": ["highs"], "std": "c99"}
    conics.append(("highs",{"highs": {"primal_feasibility_tolerance":1e-7,"solver":"choose","output_flag":False,"ipm_iteration_limit":50000}},{"quadratic": True, "dual": True, "soc": False, "codegen": codegen, "discrete": True, "sos":False}))



print(conics)


class ConicTests(casadiTestCase):

  def test_opti(self):
    for conic, qp_options, aux_options in conics:
      if conic in ["qrqp"]: continue
      opti = Opti('conic')
      x = opti.variable()
      opti.minimize(x)
      opti.subject_to(-2<=(x<=2))
      opti.solver(conic,qp_options)

      sol = opti.solve()
      self.checkarray(sol.value(x),-2,digits=5)

  def test_sos(self):

    H = DM(4,4)
    G = DM([-1,-2,-3,-1])
    A = DM([[-1,1,1,10],[1,-3,1,0],[0,1,0,-3.5]])
    LBA = DM([-inf,-inf,0])
    UBA = DM([20,30,0])
    discrete = [False,True,True,True]
    LBX = DM([0,0,0,2])
    UBX = DM([40,inf,inf,3])


    sos_groups = [[2,3]]
    sos_weights = [[25.0,18.0]]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_sos",conic,qp_options)

      options = dict(qp_options)
      options["discrete"] = discrete
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([29, 7, 22, 2]))

      if not aux_options["sos"]: continue

      options["sos_groups"] = sos_groups
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)
      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([40, 7, 0, 2]))

      options["sos_weights"] = sos_weights
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)
      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([40, 7, 0, 2]))
      print(solver_out["cost"])

      self.check_serialize(solver,solver_in)

  def test_milp(self):
    # From https://www.cs.upc.edu/~erodri/webpage/cps/theory/lp/milp/slides.pdf
    H = DM(2,2)
    G = DM([-1,-1])
    A =  DM([[-2,2],[-8,10]])

    LBA = DM([1,-inf])
    UBA = DM([inf,13])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    discrete = [True,True]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_milp",conic,qp_options)
      options = dict(qp_options)
      options["discrete"] = discrete
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([1,2]))
      self.checkarray(solver_out["cost"],DM([-3]))

      self.check_serialize(solver,solver_in)

    H = DM(3,3)
    G = DM([-1,-1,-5])
    A =  DM([[-2,2,0],[-8,10,0],[0,0,7]])

    LBA = DM([1,-inf,-4.0*7])
    UBA = DM([inf,13,4.0*7])

    LBX = DM([0]*2+[-5])
    UBX = DM([inf]*2+[5])

    discrete = [True,True,False]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_milp",conic,qp_options)

      options = dict(qp_options)
      options["discrete"] = discrete
      solver = casadi.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)
      self.checkarray(solver_out["x"],DM([1,2,4]))
      self.checkarray(solver_out["cost"],DM([-23]))

  def test_general_unconstrained(self):
    H = sparsify(DM([[1,0],[0,1]]))
    G = DM([-0.7,-2.3])
    A =  DM.zeros(0,2)

    LBA = DM.zeros(0,1)
    UBA = DM.zeros(0,1)

    LBX = DM([-inf]*2)
    UBX = DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_general_unconstrained",conic,qp_options)

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
          self.assertTrue(solver.stats()["success"])
      except:
          raise Exception(str(conic))

      self.assertAlmostEqual(solver_out["x"][0],0.7,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],2.3,max(1,6-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

      self.check_serialize(solver,solver_in)

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
      print("test_general_convex_dense",conic,qp_options)

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
      self.assertTrue(solver.stats()["success"])

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([3+1.0/9,4.0/9,0]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8-2.0/9,max(1,6-less_digits),str(conic))

      solver_in["h"]=H*4

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-6,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(conic),digits=max(1,2-less_digits))

      solver_in["h"]=0

      if 'qcqp' in str(conic): continue # Singular hessian

      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-9-1.0/3,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([10.0/3,4.0/3,0]),str(conic),digits=max(1,4-less_digits))

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

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],2,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],6,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([0,0,0]),str(conic),digits=max(1,4-less_digits))


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
      print("test_general_convex_sparse",conic,qp_options)

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

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0,-0.339076,-10.0873907,-0.252185]),6,str(conic),digits=max(1,6-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([0,2.52184767]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-6.264669320767,max(1,6-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

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
      if not("cplex" in str(conic)):
        continue
      print("test_general_nonconvex_dense",conic,qp_options)
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
      if not aux_options["quadratic"]: continue
      if "ooqp" in str(conic) or "proxqp" in str(conic):
        continue
      print("test_equality",conic,qp_options)
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

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],4.75,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([0,2,0]),str(conic),digits=max(1,6-less_digits))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
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

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,5-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,5-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([3.2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8.4,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
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
      if "qrqp" in conic: continue
      if not aux_options["quadratic"]: continue
      qp_options["dump_in"] = True
      if 'qcqp' in str(conic): continue
      print("test_degenerate_hessian",conic,qp_options)
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

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0,-2.5]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([1.5]),str(conic),digits=max(1,4-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-38.375,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

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
      print("test_no_inequality",conic,qp_options)
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

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))


      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([3.5]),str(conic),digits=max(1,6-less_digits))

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
      print("test_no_A",conic,qp_options)
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

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-34,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

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
      print("test_standard_form",conic,qp_options)
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

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([3.4]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-5.1,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

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
      print(conic,qp_options)
      if not aux_options["quadratic"]: continue
      if 'cplex' in str(conic):
        continue
      if 'worhp' in str(conic): # works but occasionaly throws segfaults, ulimit on travis?
        continue
      if 'superscs' in str(conic):
        continue
      print("test_badscaling",conic,qp_options)
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
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM.zeros(N,1),str(conic),digits=max(1,4-less_digits))

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
        print("test_redundant",conic,qp_options)
        solver = casadi.conic("qpsol",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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
        if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0,0]),str(conic),digits=max(1,6-less_digits))
        if aux_options["dual"]: self.checkarray(mtimes(A.T,solver_out["lam_a"]),DM([3.876923073076,2.4384615365384965,-1]),str(conic),digits=max(1,6-less_digits))

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
      print("test_linear",conic,qp_options)
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

      # Known bug #3038
      if conic=="superscs" and platform.system()=="Darwin":
        pass
      else:
        self.check_serialize(solver,solver_in)
      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],DM([0.5,1.5]),str(conic),digits=max(1,5-less_digits))
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,0]),str(conic),digits=max(1,5-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([0.5,-1.5,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],2.5,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])


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
      print("test_linear2",conic,qp_options)
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
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],DM([0,-3]),str(conic),digits=max(1,5-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],DM([2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],7,max(1,5-less_digits),str(conic))

  def test_overconstrained(self):
    x=SX.sym("x")
    qp={'x':x, 'f':(x-1)**2, 'g':vertcat(*[x,x,x])}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_overconstrained",conic,qp_options)
      d= dict(qp_options)
      solver = qpsol("mysolver", conic, qp, d)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10, -10, -10]
      solver_in["ubg"]=[10, 10, 10]
      solver_in["x0"]=[0]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,9,str(conic) )
      self.assertAlmostEqual(solver_out["x"][0],1,5,str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  @requires_conic("qrqp")
  def test_qrqp_prints(self):
    x = MX.sym("x")
    qp = {"x":x,"f":(x-1)**2}
    solver = qpsol("solver","qrqp",qp)
    with self.assertOutput(["last_tau","Converged"],[]):
        solver()
    if args.run_slow:
        F,_ = self.check_codegen(solver,{},std="c99",opts={})
        with self.assertOutput([],["last_tau","Converged"]): # By default, don't print
            F()
        F,_ = self.check_codegen(solver,{},std="c99",opts={"verbose_runtime":True})
        #with self.assertOutput(["last_tau","Converged"],[]): # Printing, but not captured by python stdout
        #    F()
    
    

  @requires_conic("hpipm")
  @requires_conic("qpoases")
  def test_hpipm(self):

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

    J = jacobian_old(F, 0, 0)
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
    options = {"hpipm":{"iter_max":100,"res_g_max":1e-10,"res_b_max":1e-10,"res_d_max":1e-10,"res_m_max":1e-10},"dump_in":True}
    solver = qpsol('solver', 'hpipm', prob,options)

    #solver = qpsol('solver', 'hpipm', prob,{"N":N,"nx":[2]*(N+1),"nu":[1]*N,"ng":[1]*(N+1),"tol":1e-12,"mu0":2,"max_iter":20})

    sol_ref = solver_ref(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    self.check_codegen(solver,dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg),std="c99",extralibs=["hpipm","blasfeo"])
    

    self.checkarray(sol_ref["x"], sol["x"])
    self.checkarray(sol_ref["lam_g"], sol["lam_g"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)
    self.checkarray(sol_ref["f"], sol["f"],digits=8)

    solver = nlpsol('solver', 'sqpmethod', prob,{"qpsol": "hpipm", "qpsol_options": options})
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    
    self.checkarray(sol_ref["x"], sol["x"])
    self.checkarray(sol_ref["lam_g"], sol["lam_g"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)
    self.checkarray(sol_ref["f"], sol["f"],digits=8)

  @requires_conic("hpipm")
  @requires_conic("qpoases")
  def test_hpipm_timevarying(self):
    def mat(a):
      def fl(a):
        return float(a) if len(a)>0 else 0
      return sparsify(DM([list(map(fl,i.split("\t"))) for i in a.split("\n") if len(i)>0]))
    def vec(a):
      return DM(list(map(float,a.split("\n"))))
    N = 2
    A = """
1	    0.2	  1	-1	0	  0	  0	   0	 0	0	0	0
-0.1	0.4	  0	0	 -1 	0	  0	   0	 0	0	0	0
0.3	  0.2	  0	0	  0 	-1	0	   0	 0	0	0	0
2	    0	  0.3	0	  0 	0	  0	   0	 0	0	0	0
1     1	  0.4	0	  0 	0	  0	   0	 0	0	0	0
0	    0	    1	4	  2	  1	  0.3	 -1	 0	0	0
0	    0 	  3	1	  0	  1	  0.2	 0	-1	0	0
0	    0	    1	1	  1	  1	  1	   0	 0	0	0
0	    0	    0 0	  0  	0 	0	   2	 4	0	-1
0	    0	    0	0	  0  	0	  0	   2	 3	1	0
0	    0  	  0	0	  0	  0	  0	   0	 0	0	3"""
    A = """
1	0.2	1	-1	0	0	0	0	0	0	0	0
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
    options = {"hpipm":{"iter_max":100,"res_g_max":1e-10,"res_b_max":1e-10,"res_d_max":1e-10,"res_m_max":1e-10}}
    #solver = conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},{"N":N,"nx":nx,"nu":nu,"ng":ng,"tol":1e-12,"mu0":2,"max_iter":20})
    solver = conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},options)
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

    self.check_codegen(solver,dict(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx),std="c99",extralibs=["hpipm","blasfeo"])
    
    solver = conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},options)
    sol = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"])

    self.checkarray(sol_ref["x"], sol["x"],digits=7)
    self.checkarray(sol_ref["lam_a"], sol["lam_a"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)

    # Solve again
    sol2 = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"])

    self.checkarray(sol["x"], sol2["x"])
    
    
    # extralibs=extralibs,extra_options=aux_options["codegen"]   
    print("codegen starts here")   
    self.check_codegen(solver,dict(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"]),std="c99",extralibs=["hpipm","blasfeo"])
        
  @requires_nlpsol("ipopt")
  def test_SOCP(self):

    for conic, qp_options, aux_options in conics:
      qp_options["verbose"] = True
      x = MX.sym("x")
      y = MX.sym("y")
      z = MX.sym("z")
      if not aux_options["soc"]: continue
      print("test_SOCP",conic,qp_options)

      #  min  2 x + y
      #
      #    ||  x-5 , y-7 ||_2 <= 4
      #
      #

      h = soc(vertcat(x-5,y-7),4)

      solver = casadi.qpsol("msyolver",conic,{'h':h,'x': vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      self.checkarray(res["x"],DM([5-8/sqrt(5),7-4/sqrt(5)]),conic,digits=7)
      self.checkarray(res["f"],10-16/sqrt(5)+7-4/sqrt(5),conic,digits=7)

      #  min  2 x + y
      #
      #    ||  x , y ||_2 <= ax+4
      #
      #


      a = 1.3

      h = soc(vertcat(x,y),a*x+4)


      solver = nlpsol("mysolver","ipopt",{"f":2*x+y,"x":vertcat(x,y),"g": dot(vertcat(x,y),vertcat(x,y))-(a*x+4)**2})
      res = solver(ubg=0)
      ref =  res["x"]

      solver = casadi.qpsol("msyolver",conic,{'h':h,'x': vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      xs = -(8*sqrt(5-a**2)+4*a**3-20*a)/(a**4-6*a**2+5)
      ys = 4*sqrt(5-a**2)/(a**2-5)

      self.checkarray(res["f"],2*xs+ys,conic,digits=5)
      self.checkarray(res["x"],vertcat(xs,ys),conic,digits=4 if conic in ["cplex","gurobi"] else 5)

      #  min  2 x + y
      #
      #    ||  x-5 , y-7 ||_2 <= 4
      #    ||  x/6 , y/5 ||_2 <= 1
      #

      h = diagcat(soc(vertcat(x-5,y-7),4),soc(vertcat(x/6,y/5),1))

      solver = casadi.qpsol("msyolver",conic,{'h':h,'x': vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      self.checkarray(res["x"],DM([1.655450403084473,4.805919456574478]),conic,digits=5)

      h = diagcat(soc(vertcat(-13*x+3*y+5*z-3,-12*x+12*y-6*z-2),-12*x-6*y+5*z-12),soc(vertcat(-3*x+6*y+2*z,1*x+9*y+2*z+3,-1*x-19*y+3*z-42),-3*x+6*y-10*z+27))
      solver = casadi.qpsol("msyolver",conic,{'h':h,'x': vertcat(x,y,z),"f": -2*x+1*y+5*z},qp_options)
      try:
          res = solver()
      except:
            print(solver.stats())

      self.checkarray(res["x"],DM([-5.0147928622,-5.766930599,-8.52180472]),conic,digits=4)

      self.assertTrue(solver.stats()["success"])

      # mimic a QP
      x = MX.sym("x",4)

      H = DM([[  2.834044009405148 ,  0.867080384259271 ,  0.396881849048015 ,  0.506784822363357],
         [0.867080384259271 ,  2.184677189537596 ,  0.725076945381028 ,  1.223678163433993],
         [0.396881849048015 ,  0.725076945381028,   2.838389028806589,   0.712607093594686],
         [0.506784822363357  , 1.223678163433993 ,  0.712607093594686 ,  3.340935020356804]])

      x0 = DM([1,2,3,4])
      f = -mtimes(x.T,mtimes(H,x0))

      [D,Lt,p] = ldl(H)

      F = mtimes(sqrt(diag(D)),DM.eye(4)+Lt)

      h = soc(vertcat(sqrt(2)*mtimes(F,x),1-y),1+y)

      solver = casadi.qpsol("msyolver",conic,{'x': vertcat(x,y),"f": y+f,"h":h},qp_options)
      res = solver(lbx=vertcat(-inf,-inf,-inf,-inf,1))

      self.checkarray(res["x"][:-1],x0,conic,digits=4)

  def test_no_success(self):

    x=SX.sym("x")
    y=SX.sym("y")

    f = x-y
    for conic, qp_options, aux_options in conics:
      print("test_no_success",conic,qp_options)
      opts = dict(qp_options)
      opts["error_on_fail"] = False
      solver = qpsol("solver",conic,{'x':vertcat(x,y), 'f':f,'g':vertcat(x+1,x-2)},opts)
      solver(x0=0,lbg=0,ubg=0,lbx=[-10,-10],ubx=[10,10])
      self.assertFalse(solver.stats()["success"])

      opts["error_on_fail"] = True
      solver = qpsol("solver",conic,{'x':vertcat(x,y), 'f':f,'g':vertcat(x+1,x-2)},opts)
      with self.assertInException("process"):
        solver(x0=0,lbg=0,ubg=0,lbx=[-10,-10],ubx=[10,10])

if __name__ == '__main__':
    unittest.main()
