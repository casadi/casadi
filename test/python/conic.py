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
import casadi as ca
from numpy import inf, pi
import casadi as c
import numpy
import unittest
from types import *
from helpers import args, casadiTestCase, jacobian_old, memory_heavy, requires_conic, requires_nlpsol
import os
import platform

conics = []

if "SKIP_IPOPT_TESTS" not in os.environ and ca.has_nlpsol("ipopt"):
  ipopt_options = {"fixed_variable_treatment":"relax_bounds",
                   "jac_c_constant":"yes",
                   "jac_d_constant":"yes",
                   "hessian_constant":"yes",
                   "tol":1e-12,
                   "print_level":0}
  conics.append(("nlpsol",{"nlpsol":"ipopt", "nlpsol_options.ipopt": ipopt_options},{"quadratic": True, "dual": True, "soc": False, "codegen": False, "binary":False, "discrete": False, "sos":False}))

if "SKIP_WORHP_TESTS" not in os.environ and ca.has_nlpsol("worhp"):
  worhp_options = {"TolOpti":1e-13}
  conics.append(("nlpsol",{"nlpsol":"worhp", "nlpsol_options.worhp": worhp_options},{"less_digits":1,"quadratic": True, "dual": False, "soc": False, "codegen": False, "binary":False, "discrete": False, "sos":False}))

if "SKIP_OOQP_TESTS" not in os.environ and ca.has_conic("ooqp"):
  conics.append(("ooqp",{},{"less_digits":1,"quadratic": True, "dual": True, "soc": False, "codegen": False, "binary":False, "discrete": False, "sos":False}))


if "SKIP_QPOASES_TESTS" not in os.environ and ca.has_conic("qpoases"):
  conics.append(("qpoases",{"printLevel":"low"},{"quadratic": True, "dual": True, "soc": False, "codegen": False,"binary":False, "discrete": False, "sos":False}))


if "SKIP_CPLEX_TESTS" not in os.environ and ca.has_conic("cplex"):
  conics.append(("cplex",{"cplex": {"CPX_PARAM_BARQCPEPCOMP": 1e-10,"CPX_PARAM_BAREPCOMP":1e-10}},{"quadratic": True, "dual": True, "soc": True, "codegen": False, "binary":True, "discrete": True, "sos": True}))

if "SKIP_OSQP_TESTS" not in os.environ and ca.has_conic("osqp"):
  options = ["-Wno-unused-variable","-Wno-strict-prototypes"]
  if os.name=='nt':
    options = []
  codegen = {"extralibs": ["osqp"], "extra_options": options, "std": "c99"}
  conics.append(("osqp",{"osqp":{"alpha":1,"eps_abs":1e-8,"eps_rel":1e-8}},{"quadratic": True, "dual": True, "codegen": codegen,"soc":False,"binary":False, "discrete":False}))

if "SKIP_SUPERSCS_TESTS" not in os.environ and ca.has_conic("superscs"):
  conics.append(("superscs",{"superscs": {"eps":1e-9,"do_super_scs":1, "verbose":0}},{"quadratic": True, "dual": False, "codegen": False,"soc":True,"binary":False, "discrete":False}))

# No solution for licensing on travis

if "SKIP_GUROBI_TESTS" not in os.environ and ca.has_conic("gurobi"):
  conics.append(("gurobi",{"gurobi": {"BarQCPConvTol":1e-9}},{"quadratic": True, "less_digits": True, "dual": True, "soc": True, "codegen": False,"binary":True, "discrete":True,"sos":True}))

if "SKIP_XPRESS_TESTS" not in os.environ and ca.has_conic("xpress"):
  # Note: OPTIMALITYTOL tighter than default (1e-6) can cause Xpress to
  # return infeasible solutions while reporting "optimal".  Leave at default.
  xpress_dir = os.environ.get("XPRESSDIR", "missing")
  print("xpress_dir")
  xpress_codegen = {"extralibs": ["xprs"], "std": "c99",
                    "extralibdirs": [os.path.join(xpress_dir, "lib")],
                    "extra_include": [os.path.join(xpress_dir, "include")]}
  conics.append(("xpress",{"xpress": {"OUTPUTLOG":0,"FEASTOL":1e-9,"BARGAPSTOP":1e-10,"BARDUALSTOP":1e-10}},{"quadratic": True, "dual": True, "soc": True, "codegen": xpress_codegen, "binary":True, "discrete": True, "sos":True}))

if "SKIP_MOSEK_TESTS" not in os.environ and ca.has_conic("mosek"):
  mosek_dir = os.environ.get("MOSEKDIR", os.path.expanduser("missing"))
  mosek_codegen = {"extralibs": ["mosek64"], "std": "c99",
                   "extralibdirs": [os.path.join(mosek_dir, "bin")],
                   "extra_include": [os.path.join(mosek_dir, "h")]}
  conics.append(("mosek",{"mosek": {"MSK_IPAR_LOG":0,
                                    "MSK_DPAR_INTPNT_QO_TOL_REL_GAP":1e-10,
                                    "MSK_DPAR_INTPNT_QO_TOL_PFEAS":1e-10,
                                    "MSK_DPAR_INTPNT_QO_TOL_DFEAS":1e-10,
                                    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP":1e-10,
                                    "MSK_DPAR_INTPNT_CO_TOL_PFEAS":1e-10,
                                    "MSK_DPAR_INTPNT_CO_TOL_DFEAS":1e-10}},{"quadratic": True, "less_digits": 1, "dual": True, "soc": True, "codegen": mosek_codegen, "binary":True, "discrete": True, "sos":False}))

# if has_conic("sqic"):
#   conics.append(("sqic",{},{}))

if "SKIP_CLP_TESTS" not in os.environ and ca.has_conic("clp"):
  conics.append(("clp",{"verbose":True},{"quadratic": False, "dual": True, "soc": False, "codegen": False, "binary":False, "discrete": False, "sos":False}))

if "SKIP_CBC_TESTS" not in os.environ and ca.has_conic("cbc"):
  conics.append(("cbc",{"verbose":True},{"quadratic": False, "dual": True, "soc": False, "codegen": False, "binary":True, "discrete": True, "sos":True}))

if "SKIP_QRQP_TESTS" not in os.environ and ca.has_conic("qrqp"):
  codegen = {"std":"c99"}
  conics.append(("qrqp",{"max_iter":20,"print_header":False,"print_iter":False},{"quadratic": True, "dual": True, "soc": False, "codegen": codegen, "binary":False, "discrete": False, "sos":False}))

if "SKIP_PROXQP_TESTS" not in os.environ and ca.has_conic("proxqp"):
  conics.append(("proxqp",{"proxqp":{"eps_abs":1e-11,"max_iter":1e4, "backend": "sparse"}}, {"quadratic": True, "dual": True, "soc": False, "codegen": False,"binary":False, "discrete":False,"sos":False}))
  conics.append(("proxqp",{"proxqp":{"eps_abs":1e-11,"max_iter":1e4, "backend": "dense"}}, {"quadratic": True, "dual": True, "soc": False, "codegen": False,"binary":False, "discrete":False,"sos":False}))

if "SKIP_QPALM_TESTS" not in os.environ and ca.has_conic("qpalm"):
  eps = 1e-8
  conics.append(("qpalm",{"qpalm":{"eps_abs":eps,"eps_rel":eps,"eps_abs_in":eps,"eps_prim_inf":eps}},{"quadratic": True, "dual": True, "soc": False, "codegen": False, "binary":False, "discrete": False, "sos":False}))

if "SKIP_HIGHS_TESTS" not in os.environ and ca.has_conic("highs"):
    highs_options = ["-Wno-unused-function"]
    if os.name=='nt':
      highs_options = []
    codegen = {"extralibs": ["highs"], "extra_options": highs_options, "std": "c99"}
    conics.append(("highs",{"highs": {"primal_feasibility_tolerance":1e-7,"solver":"choose","output_flag":False,"ipm_iteration_limit":50000}},{"quadratic": True, "dual": True, "soc": False, "codegen": codegen, "binary":False, "discrete": True, "sos":False}))

if "SKIP_DAQP_TESTS" not in os.environ and ca.has_conic("daqp"):
    codegen = {"extralibs": ["daqp"], "std": "c99"}
    conics.append(("daqp",{"daqp": {}},{"quadratic": True, "dual": True, "soc": False, "codegen": codegen, "binary":True, "discrete": False, "sos":False}))

if "SKIP_PIQP_TESTS" not in os.environ and ca.has_conic("piqp"):
  conics.append(("piqp",{"piqp":{"kkt_solver": "sparse_ldlt"}}, {"quadratic": True, "dual": False, "soc": False, "codegen": False,"binary":False,"discrete":False,"sos":False}))
  conics.append(("piqp",{"piqp":{"kkt_solver": "dense_cholesky"}}, {"quadratic": True, "dual": False, "soc": False, "codegen": False,"binary":False,"discrete":False,"sos":False}))


fatrop_flags = []
if os.name != 'nt':
    fatrop_flags = ["-Wno-strict-prototypes"]

print(conics)


class ConicTests(casadiTestCase):

  def test_daqp_warm_start(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    # A binary problem whose continuous relaxation is fractional.
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    inputs = dict(h=ca.DM.eye(2), g=[-0.8, -0.3], a=ca.DM([[1, 1]]),
                  lbx=[0, 0], ubx=[1, 1], lba=1, uba=1)
    for discrete in ([False, False], [True, True]):
      cold = ca.conic("cold", "daqp", structure,
                      {"discrete": discrete})
      expected = cold(**inputs)
      warm = ca.conic("warm", "daqp", structure,
                      {"discrete": discrete, "warm_start": True,
                       "daqp": {"time_limit": 0}})
      # Cover primal, dual, combined and infeasible incumbent starts.
      for start in ({"x0": expected["x"]},
                    {"lam_x0": expected["lam_x"], "lam_a0": expected["lam_a"]},
                    {"x0": expected["x"], "lam_x0": expected["lam_x"],
                     "lam_a0": expected["lam_a"]},
                    {"x0": ca.DM([0.5, 0.5])}):
        arguments = dict(inputs, **start)
        actual = warm(**arguments)
        self.checkarray(actual["x"], expected["x"])
        self.checkarray(actual["cost"], expected["cost"])
        self.assertTrue(warm.stats()["success"])
      arguments = dict(inputs, x0=expected["x"])
      self.check_serialize(warm, arguments)
      self.check_codegen(warm, arguments, extralibs=["daqp"], std="c99")

  def test_daqp_workspace_reuse(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    initial = dict(h=ca.DM.eye(2), g=[-0.8, -0.3], a=ca.DM([[1, 1]]),
                   lbx=[0, 0], ubx=[1, 1], lba=0.9, uba=1.1)
    for binary in (False, True):
      for reduction in (0, 1):
        options = {"discrete": [binary, binary], "daqp": {"eq_reduction": reduction}}
        reused = ca.conic("reused", "daqp", structure,
                          options)
        # Reallocated input objects and changing work buffers must not invalidate the cache.
        changes = ({}, {}, {"g": [-1.2, -0.4]}, {"lba": 1, "uba": 1},
                   {"a": ca.DM([[1, 0.5]]), "lba": 0.4, "uba": 1.2},
                   {"h": ca.DM([[2, 0.2], [0.2, 1]])},
                   {"h": ca.DM.eye(2)}, {"g": [-0.2, -0.2], "lba": 0, "uba": 2},
                   {"g": [-1.2, -0.4], "lba": 1, "uba": 1},
                   {"lba": 1-1e-12, "uba": 1}, {"lba": 0, "uba": 2})
        codegen_solver = ca.conic("reused_codegen", "daqp", structure,
                                 options)
        generated = self.check_codegen(codegen_solver, initial, extralibs=["daqp"], std="c99")
        if not binary:
          changes += ({"h": ca.DM([[1, 0], [0, 0]])},
                      {"h": ca.DM.eye(2)})
        arguments = dict(initial)
        for i, change in enumerate(changes):
          arguments.update(change)
          cold = ca.conic("cold", "daqp", structure, options)
          expected = cold(**arguments)
          actual = reused(**arguments)
          self.checkarray(actual["x"], expected["x"])
          self.checkarray(actual["cost"], expected["cost"])
          self.assertTrue(reused.stats()["success"])
          self.assertEqual(reused.stats()["workspace_reused"], i > 0)
          if i in (1, 2, 3):
            # No Hessian factorization or constraint-matrix update for these calls.
            self.assertEqual(reused.stats()["daqp_update_mask"] & 3, 0)
          if generated:
            result = generated["F"](**arguments)
            self.checkarray(result["x"], expected["x"])
            self.checkarray(result["cost"], expected["cost"])
        self.check_serialize(reused, arguments)
        self.check_thread_safety(reused, arguments, N=4)

  def test_daqp_reuse_equality_hessian_changes(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(2, 2)}
    solver = ca.conic("reuse_equality", "daqp", structure,
                      {"error_on_fail": False,
                       "daqp": {"eq_reduction": 1}})
    inputs = dict(h=ca.DM.eye(2), g=[-0.8, -0.3], a=ca.DM([[1, 1], [0, 0]]),
                  lbx=[0, 0], ubx=[1, 1], lba=[1, 0], uba=[1, 0])
    generated = self.check_codegen(solver, inputs, extralibs=["daqp"], std="c99")
    evaluators = [solver] + ([generated["F"]] if generated else [])
    for evaluate in evaluators:
      for unused in range(3):
        # Equality reduction must refresh saved pointers when H switches representation.
        self.checkarray(evaluate(**inputs)["x"], ca.DM([0.75, 0.25]))
        dense = dict(inputs, h=ca.DM([[4, 1], [1, 2]]), g=[-2, 0.5])
        result = evaluate(**dense)
        self.checkarray(result["x"], ca.DM([0.875, 0.125]))
        self.checkarray(result["cost"], ca.DM([-0.03125]))
      impossible = dict(inputs, lba=[1, 1], uba=[1, 1])
      if evaluate is solver:
        evaluate(**impossible)
        self.assertFalse(solver.stats()["success"])
        # DAQP >= 0.11 reports an inconsistent dependent equality as
        # overdetermined (-6) when found while forming the working set.
        self.assertIn(solver.stats()["return_status"], (-1, -6))
      else:
        with self.assertRaises(RuntimeError):
          evaluate(**impossible)
      # A failed solve must not poison the next workspace.
      self.checkarray(evaluate(**inputs)["x"], ca.DM([0.75, 0.25]))

  def test_daqp_reuse_unrelated_problems(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    rng = numpy.random.default_rng(48271)
    n, m = 8, 5
    structure = {"h": ca.Sparsity.dense(n, n), "a": ca.Sparsity.dense(m, n)}
    for reduction in (-1, 0, 1):
      for binary in (False, True):
        nb = 3 if binary else 0
        options = {"discrete": [True]*nb+[False]*(n-nb),
                   "daqp": {"eq_reduction": reduction}}
        reused = ca.conic("reuse_random", "daqp", structure,
                          options)
        for k in range(30):
          b = rng.normal(size=(n, n))
          h = numpy.diag(rng.uniform(0.2, 3, n)) if k%2 else b.T@b+numpy.eye(n)
          a = rng.normal(size=(m, n))
          witness = rng.uniform(-0.7, 0.7, n)
          witness[:nb] = rng.integers(0, 2, nb)
          center = a@witness
          lower, upper = center-0.4, center+0.6
          if k%3 == 0:
            lower[:2] = upper[:2] = center[:2]
          inputs = dict(h=h, a=a, g=-h@rng.normal(size=n),
                        lbx=[0]*nb+[-1]*(n-nb), ubx=[1]*n, lba=lower, uba=upper)
          cold = ca.conic("cold", "daqp", structure, options)
          expected, actual = cold(**inputs), reused(**inputs)
          self.assertTrue(reused.stats()["success"])
          self.checkarray(actual["x"], expected["x"])
          self.checkarray(actual["cost"], expected["cost"])
          x = numpy.array(actual["x"]).ravel()
          self.assertLessEqual(float(numpy.max(lower-a@x)), 1e-6)
          self.assertLessEqual(float(numpy.max(a@x-upper)), 1e-6)

  def test_daqp_reuse_limit_recovery(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    rng = numpy.random.default_rng(7935)
    n, m = 64, 40
    b = rng.normal(size=(n, n))/numpy.sqrt(n)
    h = b.T@b+numpy.eye(n)
    inputs = dict(h=h, a=rng.normal(size=(m, n))/numpy.sqrt(n),
                  g=-h@rng.normal(size=n)*2, lbx=[-1]*n, ubx=[1]*n,
                  lba=[-0.5]*m, uba=[0.5]*m)
    for settings, status in (({"iter_limit": 3}, -4), ({"time_limit": 1e-9}, -7)):
      solver = ca.conic("limited", "daqp",
                        {"h": ca.Sparsity.dense(n, n), "a": ca.Sparsity.dense(m, n)},
                        {"error_on_fail": False, "daqp": settings})
      solver(**inputs)
      self.assertEqual(solver.stats()["return_status"], status)
      self.assertEqual(solver.stats()["unified_return_status"], "SOLVER_RET_LIMITED")
      easy = dict(inputs, h=numpy.eye(n), g=numpy.zeros(n), a=numpy.zeros((m, n)),
                  lba=[-1]*m, uba=[1]*m)
      self.checkarray(solver(**easy)["x"], ca.DM.zeros(n))
      self.assertTrue(solver.stats()["success"])
      self.assertFalse(solver.stats()["workspace_reused"])

  def test_daqp_cold_workspace_reuse(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    rng = numpy.random.default_rng(81375)
    n, m = 8, 5
    structure = {"h": ca.Sparsity.dense(n, n), "a": ca.Sparsity.dense(m, n)}
    for binary in (False, True):
      options = {"discrete": [binary]*3+[False]*(n-3), "error_on_fail": False,
                 "daqp": {"iter_limit": 10}}
      # Defaults reuse allocations but do not reuse the previous solve's state.
      reused = ca.conic("cold_reuse", "daqp", structure, options)
      b = rng.normal(size=(n, n))
      h, a = b.T@b+numpy.eye(n), rng.normal(size=(m, n))
      for k in range(20):
        if k%4 == 0:
          b = rng.normal(size=(n, n))
          h, a = b.T@b+numpy.eye(n), rng.normal(size=(m, n))
        inputs = dict(h=h, a=a, g=-h@rng.normal(size=n),
                      lbx=[0]*3+[-1]*(n-3) if binary else [-1]*n, ubx=[1]*n,
                      lba=[-0.5]*m, uba=[0.5]*m)
        fresh = ca.conic("fresh", "daqp", structure, options)
        expected, actual = fresh(**inputs), reused(**inputs)
        self.assertEqual(reused.stats()["return_status"], fresh.stats()["return_status"])
        self.assertEqual(reused.stats()["iter_count"], fresh.stats()["iter_count"])
        if fresh.stats()["success"]:
          self.checkarray(actual["x"], expected["x"])
          self.checkarray(actual["cost"], expected["cost"])
      # Codegen exercises the same cold reset on repeated evaluations.
      unlimited = ca.conic("cold_codegen", "daqp", structure,
                            {"discrete": options["discrete"]})
      generated = self.check_codegen(unlimited, inputs, extralibs=["daqp"], std="c99")
      if generated:
        for scale in (0.5, 1, 0.7, 1.2):
          changed = dict(inputs, g=inputs["g"]*scale)
          self.checkarray(generated["F"](**changed)["x"], unlimited(**changed)["x"])

  def test_daqp_zero_row_infeasibility(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(1, 1), "a": ca.Sparsity.dense(1, 1)}
    inputs = dict(h=1, g=-0.25, a=0, lbx=-1, ubx=1, lba=-1, uba=1)
    for previous in (False, True):
      for reduction in (-1, 0, 1):
        solver = ca.conic("zero_row", "daqp", structure,
                          {"error_on_fail": False, "warm_start_previous": previous,
                           "daqp": {"eq_reduction": reduction}})
        generated = self.check_codegen(solver, inputs, extralibs=["daqp"], std="c99")
        for evaluate in [solver] + ([generated["F"]] if generated else []):
          self.checkarray(evaluate(**inputs)["x"], ca.DM([0.25]))
          impossible = dict(inputs, lba=1, uba=1)
          if evaluate is solver:
            evaluate(**impossible)
            self.assertEqual(solver.stats()["return_status"], -1)
            self.assertFalse(solver.stats()["success"])
          else:
            with self.assertRaises(RuntimeError):
              evaluate(**impossible)
          self.checkarray(evaluate(**inputs)["x"], ca.DM([0.25]))

  def test_daqp_invalid_binary_bounds_recovery(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    options = {"discrete": [True, True]}
    solver = ca.conic("binary_recovery", "daqp", structure, options)
    inputs = dict(h=ca.DM.eye(2), g=[-0.8, -0.3], a=ca.DM([[1, 1]]),
                  lbx=[0, 0], ubx=[1, 1], lba=0, uba=2)
    generated = self.check_codegen(solver, inputs, extralibs=["daqp"], std="c99")
    changed = dict(inputs, h=ca.DM.eye(2)*3, g=[-1, -2], a=ca.DM([[1, 0.5]]))
    expected = ca.conic("fresh", "daqp", structure, options)(**changed)
    for evaluate in [solver] + ([generated["F"]] if generated else []):
      evaluate(**inputs)
      with self.assertRaises(RuntimeError):
        evaluate(**dict(changed, lbx=[0.1, 0]))
      # The valid call repeats the rejected call's H/A: stale cache values
      # must not suppress recomputation of the old workspace factors.
      actual = evaluate(**changed)
      self.checkarray(actual["x"], expected["x"])
      self.checkarray(actual["cost"], expected["cost"])
      if evaluate is solver:
        self.assertFalse(solver.stats()["workspace_reused"])

  def test_daqp_sparse_matrix_updates(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    h = ca.DM([[2, 0.2, 0], [0.2, 2, 0.3], [0, 0.3, 2]])
    a = ca.DM([[1, 0, 0.5], [0, 1, -0.2]])
    structure = {"h": ca.sparsify(h).sparsity(), "a": ca.sparsify(a).sparsity()}
    inputs = dict(h=h, a=a, g=[-1, -0.5, 0.2], lbx=[-1]*3, ubx=[1]*3,
                  lba=[-0.2]*2, uba=[0.2]*2)
    solver = ca.conic("sparse_updates", "daqp", structure, {})
    generated = self.check_codegen(solver, inputs, extralibs=["daqp"], std="c99")
    changes = ({}, {}, {"a": a*0}, {"a": a},
               {"h": ca.DM.eye(3)*2}, {"h": h}, {"a": a*0.7}, {})
    arguments = dict(inputs)
    for k, change in enumerate(changes):
      arguments.update(change)
      # New caller-owned objects must not affect value-based cache detection.
      arguments["h"], arguments["a"] = ca.DM(arguments["h"]), ca.DM(arguments["a"])
      fresh = ca.conic("fresh", "daqp", structure, {})
      expected, actual = fresh(**arguments), solver(**arguments)
      self.checkarray(actual["x"], expected["x"])
      self.checkarray(actual["cost"], expected["cost"])
      if k in (1, 7):
        self.assertEqual(solver.stats()["daqp_update_mask"] & 3, 0)
      if generated:
        self.checkarray(generated["F"](**arguments)["x"], expected["x"])

  def test_daqp_bound_classification(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    inputs = dict(h=ca.DM.eye(2), a=ca.DM([[1, 1]]), g=[-0.8, -0.3],
                  lbx=[0, 0], ubx=[1, 1], lba=0, uba=2)
    for previous in (False, True):
      for reduction in (-1, 0, 1):
        options = {"warm_start_previous": previous, "error_on_fail": False,
                   "daqp": {"eq_reduction": reduction, "zero_tol": 1e-9}}
        solver = ca.conic("bounds", "daqp", structure, options)
        generated = self.check_codegen(solver, inputs, extralibs=["daqp"], std="c99")
        # Equalities are inferred by DAQP, then removed when bounds widen.
        changes = ({}, {"lba": 1, "uba": 1+1e-10},
                   {"h": ca.DM([[2, 0.2], [0.2, 1]])},
                   {"lba": 0, "uba": 2},
                   {"lbx": [0.4, 0], "ubx": [0.4+1e-10, 1]}, {},
                   {"lbx": [0, 0], "ubx": [1, 1]},
                   {"lba": 3, "uba": 4}, {"lba": 0, "uba": 2})
        arguments = dict(inputs)
        for change in changes:
          arguments.update(change)
          fresh = ca.conic("fresh", "daqp", structure, options)
          expected, actual = fresh(**arguments), solver(**arguments)
          self.assertEqual(solver.stats()["return_status"], fresh.stats()["return_status"], (previous, reduction, change))
          if fresh.stats()["success"]:
            self.checkarray(actual["x"], expected["x"])
            self.checkarray(actual["cost"], expected["cost"])
            if generated:
              self.checkarray(generated["F"](**arguments)["x"], expected["x"])
          elif generated:
            with self.assertRaises(RuntimeError):
              generated["F"](**arguments)

  def test_daqp_previous_warm_start(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    cold = ca.conic("cold", "daqp", structure, {})
    hot = ca.conic("hot", "daqp", structure, {"warm_start_previous": True})
    inputs = dict(h=ca.DM.eye(2), a=ca.DM([[1, 1]]), g=[-0.8, -0.3],
                  lbx=[0, 0], ubx=[1, 1], lba=1.2, uba=ca.inf)
    hot(**inputs)
    expected, actual = cold(**inputs), hot(**inputs)
    self.checkarray(actual["x"], expected["x"])
    self.assertLess(hot.stats()["iter_count"], cold.stats()["iter_count"])
    generated = self.check_codegen(hot, inputs, extralibs=["daqp"], std="c99")
    if generated:
      for g in ([-0.82, -0.31], [-0.7, -0.1]):
        arguments = dict(inputs, g=g)
        self.checkarray(generated["F"](**arguments)["x"], cold(**arguments)["x"])
    self.check_serialize(hot, inputs)

  def test_daqp_cold_reuse_proximal(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    structure = {"h": ca.Sparsity.dense(2, 2), "a": ca.Sparsity.dense(1, 2)}
    reused = ca.conic("reused", "daqp", structure, {})
    # The null direction makes the proximal initial center affect iteration count.
    inputs = dict(h=ca.DM([[1, 0], [0, 0]]), a=ca.DM([[1, 1]]),
                  lbx=[-1, -1], ubx=[1, 1], lba=-1, uba=1)
    for g in ([0.5, 1], [0.5, -1], [0.5, -1], [-0.5, 1]):
      arguments = dict(inputs, g=g)
      fresh = ca.conic("fresh", "daqp", structure, {})
      expected, actual = fresh(**arguments), reused(**arguments)
      self.checkarray(actual["x"], expected["x"])
      self.checkarray(actual["cost"], expected["cost"])
      self.assertEqual(reused.stats()["iter_count"], fresh.stats()["iter_count"])

  def test_daqp_setup_failure(self):
    if "SKIP_DAQP_TESTS" in os.environ or not ca.has_conic("daqp"):
      return
    solver = ca.conic("solver", "daqp", {"h": ca.Sparsity.dense(1, 1)},
                      {"error_on_fail": False, "daqp": {"eps_prox": 0}})
    solver(h=1, g=-0.5, lbx=-1, ubx=1)
    for unused in range(2):
      solver(h=-1, g=0, lbx=-1, ubx=1)
      self.assertFalse(solver.stats()["success"])
      self.assertEqual(solver.stats()["return_status"], -5)
    result = solver(h=1, g=-0.5, lbx=-1, ubx=1)
    self.assertTrue(solver.stats()["success"])
    self.checkarray(result["x"], ca.DM([0.5]))

  @memory_heavy()
  def test_binary(self):
    ca.DM.rng(1)
    n = 10
    m = 8
    ms = 5
    x = ca.MX.sym("x", n)
    # --- Hessian: dense symmetric matrix ---
    M = ca.DM.rand(n,n)
    H = 0.5 * (M @ M.T)  # symmetrize

    # --- Linear term ---
    f = 100 * ca.DM.rand(n,1)

    # --- Constraint matrix and bounds ---
    A = ca.DM.rand(m,n)
    bupper = 20 * ca.DM.rand(m,1)
    blower = -20 * ca.DM.rand(m,1)

    solver_in = dict(lbx=[0] * ms + [-10] * (n-ms), ubx=[1] * ms + [10] * (n-ms),
                 lbg=blower, ubg=bupper)
    for conic, qp_options, aux_options in conics:
      if not aux_options["binary"]: continue
      
      if conic not in ["daqp"]:
        solver = ca.qpsol('solver', conic,
                          {'f': f.T @ x, 'x': x, "g": A @ x},
                          {'discrete': [1] * ms + [0] * (n-ms)}
                          )
                          
        solver(**solver_in)

        if aux_options["codegen"]:
          self.check_codegen(solver,solver.convert_in(solver_in),**aux_options["codegen"])        
                   
      if not aux_options["quadratic"]: continue
      if conic not in ["highs"]:

        # Using casadi qpsol
        solver = ca.qpsol('solver', conic,
                          {'f': 0.5*(x.T @ H @ x) + f.T @ x, 'x': x, "g": A @ x},
                          {'discrete': [1] * ms + [0] * (n-ms)}
                          )
        solver(**solver_in)

        digits = 15
        if conic == "mosek" and platform.system()=="Darwin":
            digits = 5
        if aux_options["codegen"]:
          self.check_codegen(solver,solver.convert_in(solver_in),digits=digits,**aux_options["codegen"])

  @requires_conic("mosek")
  def test_mosek_codegen_option_precision(self):
    import re
    options = {"MSK_DPAR_MIO_TOL_ABS_RELAX_INT": 1.234567890123456e-7,
               "MSK_DPAR_OPTIMIZER_MAX_TIME": 1234567890.1234567}
    solver = ca.conic("solver", "mosek", {"h": ca.Sparsity.dense(1, 1)},
                      {"mosek": options})
    cg = ca.CodeGenerator("mosek_precision.c")
    cg.add(solver)
    emitted = dict(re.findall(r'MSK_putnadouparam\(d->task, "([^"]+)", ([^)]+)\)',
                              cg.dump()))
    self.assertEqual(set(emitted), set(options))
    for name, value in options.items():
      self.assertEqual(float(emitted[name]), value, name)

  def test_opti(self):
    for conic, qp_options, aux_options in conics:
      if conic in ["qrqp"]: continue
      if conic in ["daqp"]: continue
      opti = ca.Opti('conic')
      x = opti.variable()
      opti.minimize(x)
      opti.subject_to(-2<=(x<=2))
      opti.solver(conic,qp_options)

      sol = opti.solve()
      self.checkarray(sol.value(x),-2,digits=5)

  def test_sos(self):

    H = ca.DM(4,4)
    G = ca.DM([-1,-2,-3,-1])
    A = ca.DM([[-1,1,1,10],[1,-3,1,0],[0,1,0,-3.5]])
    LBA = ca.DM([-inf,-inf,0])
    UBA = ca.DM([20,30,0])
    discrete = [False,True,True,True]
    LBX = ca.DM([0,0,0,2])
    UBX = ca.DM([40,inf,inf,3])


    sos_groups = [[2,3]]
    sos_weights = [[25.0,18.0]]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_sos",conic,qp_options)

      options = dict(qp_options)
      options["discrete"] = discrete
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],ca.DM([29, 7, 22, 2]))

      if not aux_options["sos"]: continue

      options["sos_groups"] = sos_groups
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)
      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],ca.DM([40, 7, 0, 2]))

      options["sos_weights"] = sos_weights
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)
      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],ca.DM([40, 7, 0, 2]))
      print(solver_out["cost"])

      self.check_serialize(solver,solver_in)

  def test_milp(self):
    # From https://www.cs.upc.edu/~erodri/webpage/cps/theory/lp/milp/slides.pdf
    H = ca.DM(2,2)
    G = ca.DM([-1,-1])
    A =  ca.DM([[-2,2],[-8,10]])

    LBA = ca.DM([1,-inf])
    UBA = ca.DM([inf,13])

    LBX = ca.DM([0]*2)
    UBX = ca.DM([inf]*2)

    discrete = [True,True]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_milp",conic,qp_options)
      options = dict(qp_options)
      options["discrete"] = discrete
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["x"],ca.DM([1,2]))
      self.checkarray(solver_out["cost"],ca.DM([-3]))

      self.check_serialize(solver,solver_in)

    H = ca.DM(3,3)
    G = ca.DM([-1,-1,-5])
    A =  ca.DM([[-2,2,0],[-8,10,0],[0,0,7]])

    LBA = ca.DM([1,-inf,-4.0*7])
    UBA = ca.DM([inf,13,4.0*7])

    LBX = ca.DM([0]*2+[-5])
    UBX = ca.DM([inf]*2+[5])

    discrete = [True,True,False]

    for conic, qp_options, aux_options in conics:
      if not aux_options["discrete"]: continue
      print("test_milp",conic,qp_options)

      options = dict(qp_options)
      options["discrete"] = discrete
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      solver_out = solver(**solver_in)
      self.checkarray(solver_out["x"],ca.DM([1,2,4]))
      self.checkarray(solver_out["cost"],ca.DM([-23]))

  def test_general_unconstrained(self):
    H = ca.sparsify(ca.DM([[1,0],[0,1]]))
    G = ca.DM([-0.7,-2.3])
    A =  ca.DM.zeros(0,2)

    LBA = ca.DM.zeros(0,1)
    UBA = ca.DM.zeros(0,1)

    LBX = ca.DM([-inf]*2)
    UBX = ca.DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_general_unconstrained",conic,qp_options)

      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([[1, 1],[-1, 2],[2, 1]])

    LBA = ca.DM([-inf]*3)
    UBA = ca.DM([2, 2, 3])

    LBX = ca.DM([0]*2)
    UBX = ca.DM([inf]*2)

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if conic == "nlpsol": continue
      if not aux_options["quadratic"]: continue
      print("test_general_convex_dense",conic,qp_options)
      if conic in ["daqp"]: continue

      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)
      try:
        solver.stats()
      except:
        pass #  nosegfault

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
      self.assertEqual(solver.stats()["unified_return_status"],"SOLVER_RET_SUCCESS")

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([3+1.0/9,4.0/9,0]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8-2.0/9,max(1,6-less_digits),str(conic))

      solver_in["h"]=H*4

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1,max(1,3-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-6,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([2,0,0]),str(conic),digits=max(1,2-less_digits))

      solver_in["h"]=0

      if 'qcqp' in str(conic): continue # Singular hessian

      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],2.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-9-1.0/3,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],0,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([10.0/3,4.0/3,0]),str(conic),digits=max(1,4-less_digits))

      solver_in["lba"]=[-inf]*3 #  Upper _and_ lower
      solver_in["uba"]=[inf]*3  #  bounds infinite?

      solver_in["ubx"]=5

      if "worhp" in str(qp_options):
        with self.assertRaises(Exception):
          solver_out = solver(solver_in)  # pyright: ignore[reportCallIssue,reportArgumentType]
        continue
      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["cost"],-40,max(1,5-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],2,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],6,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([0,0,0]),str(conic),digits=max(1,4-less_digits))


  def test_general_convex_sparse(self):
    self.message("Convex sparse QP with solvers: " + str([conic for conic,options,aux_options in conics]))
    H = c.diag([2,1,0.2,0.7,1.3])

    H[1,2]=0.1
    H[2,1]=0.1

    G = ca.DM([-2,-6,1,0,0])
    A =  ca.DM([[1, 0,0.1,0.7,-1],[0.1, 2,-0.3,4,0.1]])
    A = ca.sparsify(A)

    LBA = ca.DM([-inf])
    UBA = ca.DM([2, 2])

    LBX = ca.DM([0]*5)
    UBX = ca.DM([inf]*5)


    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_general_convex_sparse",conic,qp_options)

      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      self.checkarray(solver_out["x"],ca.DM([0.873908,0.95630465,0,0,0]),str(conic),digits=max(1,6-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0,-0.339076,-10.0873907,-0.252185]),str(conic),digits=max(1,6-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([0,2.52184767]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-6.264669320767,max(1,6-less_digits),str(conic))

      if aux_options["codegen"]:
        # DAQP equality inference can change roundoff between fresh/reused calls.
        self.check_codegen(solver,solver_in,digits=13 if conic=="daqp" else 15,
                           **aux_options["codegen"])

  def test_general_nonconvex_dense(self):
    self.message("Non convex dense QP with solvers: " + str([conic for conic,options,aux_options in conics]))
    H = ca.DM([[1,-1],[-1,-2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([[1, 1],[-1, 2],[2, 1]])

    LBA = ca.DM([-inf]*3)
    UBA = ca.DM([2, 2, 3])

    LBX = ca.DM([0]*2)
    UBX = ca.DM([inf]*2)

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_general_nonconvex_dense",conic,qp_options)
      solver = ca.conic("mysolver",conic, {'h':H.sparsity(),'a':A.sparsity()},qp_options)

      solver_in = {}
      solver_in["h"]=H
      solver_in["g"]=G
      solver_in["a"]=A
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lba"]=LBA
      solver_in["uba"]=UBA

      self.assertRaises(Exception,lambda : solver(solver_in))  # pyright: ignore[reportCallIssue,reportArgumentType]

  @memory_heavy()
  def test_equality(self):
    self.message("Regression 452 test: equality constraints give wrong multipliers")
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])

    options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      if "ooqp" in str(conic) or "proxqp" in str(conic):
        continue
      print("test_equality",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':ca.Sparsity.dense(3,2)},qp_options)

      try:
        less_digits=aux_options["less_digits"]
      except:
        less_digits=0

      A =  ca.DM([[1, 1],[-1, 2],[2, 1]])

      LBA = ca.DM([-inf]*3)
      UBA = ca.DM([2, 2, 3])

      LBX = ca.DM([0.5,0])
      UBX = ca.DM([0.5,inf])

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
          solver_out = solver(solver_in)  # pyright: ignore[reportCallIssue,reportArgumentType]
        continue

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,max(1,6-less_digits),str(conic))
      self.assertAlmostEqual(solver_out["x"][1],1.25,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][0],4.75,max(1,6-less_digits),str(conic))
      if aux_options["dual"]: self.assertAlmostEqual(solver_out["lam_x"][1],0,max(1,6-less_digits),str(conic))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([0,2,0]),str(conic),digits=max(1,6-less_digits))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
      self.assertAlmostEqual(solver_out["cost"][0],-7.4375,max(1,6-less_digits),str(conic))

      A =  ca.DM([[1, 1],[-1, 2],[2, 1]])
      LBA = ca.DM([2,-inf,-inf])
      UBA = ca.DM([2, inf, inf])

      LBX = ca.DM([-inf]*2)
      UBX = ca.DM([inf]*2)


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

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([3.2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-8.4,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
  @memory_heavy()
  def test_degenerate_hessian(self):
    self.message("Degenerate hessian")

    H = ca.DM([[1,-1,0],[-1,2,0],[0,0,0]])
    H = ca.sparsify(H)
    G = ca.DM([-2,-6,1])
    A =  ca.DM([[1, 1,1]])




    LBA = ca.DM([0.5])
    UBA = ca.DM([0.5])

    LBX = ca.DM([-10])
    UBX = ca.DM([10])

    for conic, qp_options, aux_options in conics:
      if "qrqp" in conic: continue
      if not aux_options["quadratic"]: continue
      qp_options["dump_in"] = True
      if 'qcqp' in str(conic): continue
      if conic in ["daqp"]: continue
      print("test_degenerate_hessian",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      self.checkarray(solver_out["x"],ca.DM([5.5,5,-10]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0,-2.5]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([1.5]),str(conic),digits=max(1,4-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-38.375,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_no_inequality(self):
    self.message("No inequalities present")
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([[1, 1]])




    LBA = ca.DM([0.5])
    UBA = ca.DM([0.5])

    LBX = ca.DM([-10])
    UBX = ca.DM([10])


    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_no_inequality",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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


      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([3.5]),str(conic),digits=max(1,6-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-3.375,max(1,6-less_digits),str(conic))

  def test_no_A(self):
    self.message("No A present")
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])
    A =  ca.DM(0,2)

    LBA = ca.DM(0,1)
    UBA = ca.DM(0,1)

    LBX = ca.DM([-10])
    UBX = ca.DM([10])



    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_no_A",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      self.checkarray(solver_out["x"],ca.DM([10,8]),str(conic),digits=max(1,3-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-34,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_standard_form(self):
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([1,1]).T

    LBA = ca.DM([-inf])
    UBA = ca.DM([1])

    LBX = ca.DM([-10])
    UBX = ca.DM([10])

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_standard_form",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)


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

      self.checkarray(solver_out["x"],ca.DM([-0.2,1.2]),str(conic),digits=max(1,3-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0]),str(conic),digits=max(1,4-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([3.4]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],-5.1,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  @memory_heavy()
  def test_badscaling(self):
    #return
    self.message("Badly scaled problem")
    N = 50
    H = c.diag(list(range(1,N+1)))
    x0 = ca.DM(list(range(N)))

    G = -1.0*(H @ x0)

    A =  ca.DM(0,N)

    LBX = ca.DM([-1000]*N)
    UBX = ca.DM([1000]*N)

    for conic, qp_options, aux_options in conics:
      print(conic,qp_options)
      if not aux_options["quadratic"]: continue
      if 'worhp' in str(conic): # works but occasionaly throws segfaults, ulimit on travis?
        continue
      if 'superscs' in str(conic):
        continue
      print("test_badscaling",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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
      self.assertAlmostEqual(solver_out["cost"][0],-0.5*ca.mtimes([x0.T,H,x0]),max(1,3-less_digits),str(conic))
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM.zeros(N,1),str(conic),digits=max(1,4-less_digits))

  def test_redundant(self):
    self.message("Redundant constraints")

    H = ca.DM([[1,-1,0],[-1,2,0],[0,0,0]])
    G = ca.DM([-2,-6,1])
    a = ca.DM([1,0,1])
    a_ = ca.DM([0,1,-2])

    for w0,w1 in [(0,2),(1,1),(0.1,0.6)]:

      A =  ca.vertcat(*[a.T,a_.T,(w0*a+w1*a_).T])

      LBA = ca.DM([0,0,0])
      UBA = ca.DM([0.5,0.3,w0*0.5+w1*0.3])

      LBX = ca.DM([-10])
      UBX = ca.DM([10])

      options = {"mutol": 1e-12, "artol": 1e-12, "tol":1e-12}

      for conic, qp_options, aux_options in conics:
        if not aux_options["quadratic"]: continue
        if 'qcqp' in str(conic): continue
        print("test_redundant",conic,qp_options)
        if conic in ["daqp"]: continue
        solver = ca.conic("qpsol",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

        self.checkarray(solver_out["x"],ca.DM([-0.19230768069,1.6846153915,0.692307690769276]),str(conic),digits=max(1,6-less_digits))
        self.assertAlmostEqual(solver_out["cost"][0],-5.850384678537,max(1,5-less_digits),str(conic))
        if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0,0]),str(conic),digits=max(1,6-less_digits))
        if aux_options["dual"]: self.checkarray(A.T @ solver_out["lam_a"],ca.DM([3.876923073076,2.4384615365384965,-1]),str(conic),digits=max(1,6-less_digits))

  def test_linear(self):
    H = ca.DM(2,2)
    A = ca.DM([ [-1,1],[1,1],[1,-2]])
    LBA = ca.DM([ -inf, 2, -inf ])
    UBA = ca.DM([ 1, inf, 4 ])
    LBX = ca.DM([ -inf, 0 ])
    UBX = ca.DM([ inf, inf ])
    G = ca.DM([ 2, 1 ])


    for conic, qp_options, aux_options in conics:
      if 'qcqp' in str(conic): continue
      if conic in ["daqp"]: continue
      print("test_linear",conic,qp_options)
      solver = ca.conic("mysolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      self.checkarray(solver_out["x"],ca.DM([0.5,1.5]),str(conic),digits=max(1,5-less_digits))
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,0]),str(conic),digits=max(1,5-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([0.5,-1.5,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],2.5,max(1,5-less_digits),str(conic))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])


  def test_linear2(self):
    H = ca.DM(2,2)
    A = ca.DM([[-1,1],[1,1],[1,-2]])
    LBA = ca.DM([ -inf, 2, -inf ])
    UBA = ca.DM([ 1, inf, 4 ])
    LBX = ca.DM([ -inf, 3 ])
    UBX = ca.DM([ inf, 3 ])
    G = ca.DM([ 2.0, 1.0 ])


    for conic, qp_options, aux_options in conics:
      if 'qcqp' in str(conic): continue
      if 'nlp' in str(conic): continue
      if conic in ["daqp"]: continue
      print("test_linear2",conic,qp_options)
      solver = ca.conic("msyolver",conic,{'h':H.sparsity(),'a':A.sparsity()},qp_options)

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

      self.checkarray(solver_out["x"],ca.DM([2,3]),str(conic),digits=max(1,5-less_digits))
      if aux_options["dual"]: self.checkarray(solver_out["lam_x"],ca.DM([0,-3]),str(conic),digits=max(1,5-less_digits))

      if aux_options["dual"]: self.checkarray(solver_out["lam_a"],ca.DM([2,0,0]),str(conic),digits=max(1,5-less_digits))

      self.assertAlmostEqual(solver_out["cost"][0],7,max(1,5-less_digits),str(conic))

  def test_overconstrained(self):
    x=ca.SX.sym("x")
    qp={'x':x, 'f':(x-1)**2, 'g':ca.vertcat(*[x,x,x])}

    for conic, qp_options, aux_options in conics:
      if not aux_options["quadratic"]: continue
      print("test_overconstrained",conic,qp_options)
      d= dict(qp_options)
      solver = ca.qpsol("mysolver", conic, qp, d)
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
    x = ca.MX.sym("x")
    qp = {"x":x,"f":(x-1)**2}
    solver = ca.qpsol("solver","qrqp",qp)
    with self.assertOutput(["last_tau","Converged"],[]):
        solver()
    if args.run_slow:
        cg = self.check_codegen(solver,{},std="c99",opts={})
        F = cg["F"]
        with self.assertOutput([],["last_tau","Converged"]): # By default, don't print
            F()
        cg = self.check_codegen(solver,{},std="c99",opts={"verbose_runtime":True})
        F = cg["F"]
        #with self.assertOutput(["last_tau","Converged"],[]): # Printing, but not captured by python stdout
        #    F()
    
    

  @requires_conic("hpipm")
  @requires_conic("qpoases")
  def test_hpipm(self):

    inf = 100
    T = 10. # Time horizon
    N = 4 # number of control intervals

    # Declare model variables
    x1 = ca.MX.sym('x1')
    x2 = ca.MX.sym('x2')
    x = ca.vertcat(x1, x2)
    u = ca.MX.sym('u')

    # Model equations
    xdot = ca.vertcat(0.6*x1 - 1.11*x2 + 0.3*u-0.03, 0.7*x1+0.01)

    # Objective term
    L = x1**2 + 3*x2**2 + 7*u**2 -0.4*x1*x2-0.3*x1*u+u -x1-2*x2

    # Fixed step Runge-Kutta 4 integrator
    F = ca.Function('F', [x, u], [x+xdot, L])

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

    Xs = ca.SX.sym('X', 2, 1, N+1)
    Us = ca.SX.sym('U', 1, 1, N+1)

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

    J+= Xs[-1].T @ Xs[-1]

    w += [Xs[-1]]
    lbw += [-inf, -inf]
    ubw += [  inf,  inf]
    w0  += [0, 1]

    # Create an NLP solver
    prob = {'f': J, 'x': ca.vertcat(*w), 'g': ca.vertcat(*g)}


    J = ca.Function("J",[prob["x"]],[ca.jacobian(prob["g"],prob["x"])])
    J(w0).print_dense()


    solver_ref = ca.qpsol('solver', 'qpoases', prob)
    options = {"hpipm":{"iter_max":100,"res_g_max":1e-10,"res_b_max":1e-10,"res_d_max":1e-10,"res_m_max":1e-10},"dump_in":True}
    solver = ca.qpsol('solver', 'hpipm', prob,options)

    #solver = qpsol('solver', 'hpipm', prob,{"N":N,"nx":[2]*(N+1),"nu":[1]*N,"ng":[1]*(N+1),"tol":1e-12,"mu0":2,"max_iter":20})

    sol_ref = solver_ref(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)

    self.check_codegen(solver,dict(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg),std="c99",extralibs=["hpipm","blasfeo"],extra_options=fatrop_flags,digits=10)
    

    self.checkarray(sol_ref["x"], sol["x"])
    self.checkarray(sol_ref["lam_g"], sol["lam_g"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)
    self.checkarray(sol_ref["f"], sol["f"],digits=8)

    solver = ca.nlpsol('solver', 'sqpmethod', prob,{"qpsol": "hpipm", "qpsol_options": options})
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
      return ca.sparsify(ca.DM([list(map(fl,i.split("\t"))) for i in a.split("\n") if len(i)>0]))
    def vec(a):
      return ca.DM(list(map(float,a.split("\n"))))
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
    solver = ca.conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},options)
    solver_ref = ca.conic('solver', 'qpoases', {"a": A.sparsity(), "h": H.sparsity()})

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

    self.check_codegen(solver,dict(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx),std="c99",extralibs=["hpipm","blasfeo"],extra_options=fatrop_flags,digits=9)
    
    solver = ca.conic('solver', 'hpipm', {"a": A.sparsity(), "h": H.sparsity()},options)
    sol = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"])

    self.checkarray(sol_ref["x"], sol["x"],digits=7)
    self.checkarray(sol_ref["lam_a"], sol["lam_a"],digits=8)
    self.checkarray(sol_ref["lam_x"], sol["lam_x"],digits=8)

    # Solve again
    sol2 = solver(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"])

    self.checkarray(sol["x"], sol2["x"])
    
    
    # extralibs=extralibs,extra_options=aux_options["codegen"]   
    print("codegen starts here")   
    self.check_codegen(solver,dict(a=A,h=H,lba=lbg,uba=ubg,g=g,lbx=lbx,ubx=ubx,x0=sol["x"],lam_a0=sol["lam_a"],lam_x0=sol["lam_x"]),std="c99",extralibs=["hpipm","blasfeo"],extra_options=fatrop_flags,digits=9)
        
  @requires_nlpsol("ipopt")
  def test_SOCP(self):

    for conic, qp_options, aux_options in conics:
      qp_options["verbose"] = True
      x = ca.MX.sym("x")
      y = ca.MX.sym("y")
      z = ca.MX.sym("z")
      if not aux_options["soc"]: continue
      print("test_SOCP",conic,qp_options)

      #  min  2 x + y
      #
      #    ||  x-5 , y-7 ||_2 <= 4
      #
      #

      h = ca.soc(ca.vertcat(x-5,y-7),4)

      solver = ca.qpsol("msyolver",conic,{'h':h,'x': ca.vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      self.checkarray(res["x"],ca.DM([5-8/ca.sqrt(5),7-4/ca.sqrt(5)]),conic,digits=7)
      self.checkarray(res["f"],10-16/ca.sqrt(5)+7-4/ca.sqrt(5),conic,digits=7)

      #  min  2 x + y
      #
      #    ||  x , y ||_2 <= ax+4
      #
      #


      a = 1.3

      h = ca.soc(ca.vertcat(x,y),a*x+4)


      solver = ca.nlpsol("mysolver","ipopt",{"f":2*x+y,"x":ca.vertcat(x,y),"g": ca.dot(ca.vertcat(x,y),ca.vertcat(x,y))-(a*x+4)**2})
      res = solver(ubg=0)
      ref =  res["x"]

      solver = ca.qpsol("msyolver",conic,{'h':h,'x': ca.vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      xs = -(8*ca.sqrt(5-a**2)+4*a**3-20*a)/(a**4-6*a**2+5)
      ys = 4*ca.sqrt(5-a**2)/(a**2-5)

      self.checkarray(res["f"],2*xs+ys,conic,digits=5)
      self.checkarray(res["x"],ca.vertcat(xs,ys),conic,digits=4 if conic in ["cplex","gurobi","xpress"] else 5)

      #  min  2 x + y
      #
      #    ||  x-5 , y-7 ||_2 <= 4
      #    ||  x/6 , y/5 ||_2 <= 1
      #

      h = ca.diagcat(ca.soc(ca.vertcat(x-5,y-7),4),ca.soc(ca.vertcat(x/6,y/5),1))

      solver = ca.qpsol("msyolver",conic,{'h':h,'x': ca.vertcat(x,y),"f": 2*x+y},qp_options)

      res = solver()

      self.checkarray(res["x"],ca.DM([1.655450403084473,4.805919456574478]),conic,digits=5)

      h = ca.diagcat(ca.soc(ca.vertcat(-13*x+3*y+5*z-3,-12*x+12*y-6*z-2),-12*x-6*y+5*z-12),ca.soc(ca.vertcat(-3*x+6*y+2*z,1*x+9*y+2*z+3,-1*x-19*y+3*z-42),-3*x+6*y-10*z+27))
      solver = ca.qpsol("msyolver",conic,{'h':h,'x': ca.vertcat(x,y,z),"f": -2*x+1*y+5*z},qp_options)
      try:
          res = solver()
      except:
            print(solver.stats())

      self.checkarray(res["x"],ca.DM([-5.0147928622,-5.766930599,-8.52180472]),conic,digits=4)

      self.assertTrue(solver.stats()["success"])

      # mimic a QP
      x = ca.MX.sym("x",4)

      H = ca.DM([[  2.834044009405148 ,  0.867080384259271 ,  0.396881849048015 ,  0.506784822363357],
         [0.867080384259271 ,  2.184677189537596 ,  0.725076945381028 ,  1.223678163433993],
         [0.396881849048015 ,  0.725076945381028,   2.838389028806589,   0.712607093594686],
         [0.506784822363357  , 1.223678163433993 ,  0.712607093594686 ,  3.340935020356804]])

      x0 = ca.DM([1,2,3,4])
      f = -(x.T @ H @ x0)

      [D,Lt,p] = ca.ldl(H)

      F = ca.sqrt(ca.diag(D)) @ (ca.DM.eye(4)+Lt)

      h = ca.soc(ca.vertcat(ca.sqrt(2)*(F @ x),1-y),1+y)

      solver = ca.qpsol("msyolver",conic,{'x': ca.vertcat(x,y),"f": y+f,"h":h},qp_options)
      res = solver(lbx=ca.vertcat(-inf,-inf,-inf,-inf,1))

      self.checkarray(res["x"][:-1],x0,conic,digits=4)

  def test_no_success(self):

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    f = x-y
    for conic, qp_options, aux_options in conics:
      print("test_no_success",conic,qp_options)
      opts = dict(qp_options)
      opts["error_on_fail"] = False
      solver = ca.qpsol("solver",conic,{'x':ca.vertcat(x,y), 'f':f,'g':ca.vertcat(x+1,x-2)},opts)
      solver(x0=0,lbg=0,ubg=0,lbx=[-10,-10],ubx=[10,10])
      self.assertFalse(solver.stats()["success"])

      opts["error_on_fail"] = True
      solver = ca.qpsol("solver",conic,{'x':ca.vertcat(x,y), 'f':f,'g':ca.vertcat(x+1,x-2)},opts)
      with self.assertInException("process"):
        solver(x0=0,lbg=0,ubg=0,lbx=[-10,-10],ubx=[10,10])

  @requires_conic("osqp")
  def test_expand(self):

    x=ca.MX.sym("x")

    qp = {'x':x, 'f':(x-1)}
    solver = ca.qpsol("solver","osqp",qp)

    solver_qp = solver.find_function("solver_qp")
    self.assertTrue(solver_qp.is_a("MXFunction"))
    
    solver = ca.qpsol("solver","osqp",qp,{"expand":True})
    
    solver_qp = solver.find_function("solver_qp")
    self.assertTrue(solver_qp.is_a("SXFunction"))
    
  @requires_conic("osqp")
  def test_postpone_expand(self):

    x=ca.MX.sym("x")
    
    J = ca.Function("jac_f", [x, ca.MX(1,1)], [0.5*x], ['x', 'out_r'], ['jac_r_x'],{"always_inline":True})
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    
    qp = {'x':x, 'f':f(x-1)}
    solver = ca.qpsol("solver","osqp",qp,{"print_time":True})
    solver(x0=2)
    solver_qp = solver.find_function("solver_qp")
    self.checkarray(solver_qp(x=2)["H"],0.5)
    self.assertTrue(solver_qp.is_a("MXFunction"))
    
    solver = ca.qpsol("solver","osqp",qp,{"expand":True,"print_time":True})
    solver(x0=2)
    solver_qp = solver.find_function("solver_qp")
    self.checkarray(solver_qp(x=2)["H"],2)
    self.assertTrue(solver_qp.is_a("SXFunction"))
    
    solver = ca.qpsol("solver","osqp",qp,{"print_time":True, "postpone_expand":True})
    solver(x0=2)
    solver_qp = solver.find_function("solver_qp")
    self.checkarray(solver_qp(x=2)["H"],0.5)
    self.assertTrue(solver_qp.is_a("MXFunction"))
    
    solver = ca.qpsol("solver","osqp",qp,{"expand":True,"print_time":True, "postpone_expand":True})
    solver(x0=2)
    solver_qp = solver.find_function("solver_qp")
    self.checkarray(solver_qp(x=2)["H"],0.5)
    self.assertTrue(solver_qp.is_a("SXFunction"))
    
        
@requires_conic("xpress")
class XpressSpecificTests(casadiTestCase):
  """Xpress-specific tests for features not covered by the generic loop."""

  def test_mip_warmstart(self):
    """mip_start=True: warm-start hint accepted, 'User solution' in log."""
    import tempfile, os

    x = ca.MX.sym("x", 3)
    # Binary MILP: minimise -sum(x) s.t. x[0]+x[1]+x[2] <= 2
    # Optimal: any two of the three = 1, cost = -2
    logfile = tempfile.mktemp(suffix=".log")
    try:
      solver = ca.qpsol("s", "xpress",
                     {"x": x, "f": -ca.sum1(x), "g": x[0] + x[1] + x[2]},
                     {"discrete": [1, 1, 1],
                      "mip_start": True,
                      "log_file": logfile,
                      "xpress": {"OUTPUTLOG": 1}})

      sol = solver(x0=ca.DM([1, 1, 0]), lbx=ca.DM([0, 0, 0]), ubx=ca.DM([1, 1, 1]),
                   lbg=-inf, ubg=2)

      self.assertTrue(solver.stats()["success"], "MIP with mip_start should succeed")
      self.assertAlmostEqual(float(sol["f"]), -2.0, places=6,
                             msg="MILP optimal cost should be -2")

      with open(logfile) as f:
        log = f.read()
      self.assertIn("User solution", log,
                    "Xpress log must contain 'User solution' when mip_start=True")
    finally:
      del solver  # release Xpress log file handle before unlinking
      if os.path.exists(logfile):
        os.unlink(logfile)


  # XPRSiisfirst has an off-by-one heap overflow (libxprs) that corrupts the
  # heap on macOS; compute_iis is opt-in and this test is skipped there.
  @unittest.skipIf(platform.system()=="Darwin", "Xpress XPRSiisfirst corrupts the heap on macOS")
  def test_iis_infeasible_lp(self):
    """IIS rows/cols are returned in get_stats() for an infeasible LP."""
    x = ca.MX.sym("x")
    # Infeasible via constraints: x >= 5 AND -x >= 5 (i.e. x <= -5)
    # Bounds keep x free so CasADi's pre-check passes.
    solver = ca.qpsol("s", "xpress",
                   {"x": x, "f": x, "g": ca.vertcat(x, -x)},
                   {"error_on_fail": False,
                    "compute_iis": True,
                    "xpress": {"OUTPUTLOG": 0}})
    sol = solver(x0=0, lbx=-1e10, ubx=1e10, lbg=ca.DM([5, 5]), ubg=ca.DM([1e10, 1e10]))

    stats = solver.stats()
    self.assertFalse(stats["success"], "Solver should report failure for infeasible LP")
    self.assertIn("iis_rows", stats,
                  "get_stats() must contain 'iis_rows' for infeasible LP")
    self.assertIn("iis_cols", stats,
                  "get_stats() must contain 'iis_cols' for infeasible LP")
    self.assertGreater(len(stats["iis_rows"]), 0,
                       "IIS must identify at least one infeasible row")
    # When no solution is available, x and duals must be NaN (not garbage).
    self.assertTrue(numpy.all(numpy.isnan(numpy.array(sol["x"]).flatten())),
                    "x must be NaN when no primal solution is available")

  @unittest.skipIf(platform.system()=="Darwin", "Xpress XPRSiisfirst corrupts the heap on macOS")
  def test_iis_infeasible_milp(self):
    """IIS is returned in get_stats() for an infeasible MIP."""
    x = ca.MX.sym("x")
    # Binary variable (0/1) with constraint x >= 2 -> infeasible
    solver = ca.qpsol("s", "xpress",
                   {"x": x, "f": x, "g": x},
                   {"discrete": [1],
                    "error_on_fail": False,
                    "compute_iis": True,
                    "xpress": {"OUTPUTLOG": 0}})
    solver(x0=0, lbx=ca.DM([0]), ubx=ca.DM([1]), lbg=2, ubg=1e10)

    stats = solver.stats()
    self.assertFalse(stats["success"])
    self.assertIn("iis_rows", stats)
    self.assertIn("iis_cols", stats)
    self.assertGreater(len(stats["iis_rows"]) + len(stats["iis_cols"]), 0,
                       "IIS must contain at least one row or column")


if __name__ == '__main__':
    unittest.main()
