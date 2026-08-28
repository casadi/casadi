#
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
import sys
from numpy import inf, pi
import casadi as c
import numpy
import numpy as np
import unittest
from types import *
from helpers import *
import itertools
import copy
from casadi.tools import capture_stdout

import os
#GlobalOptions.setCatchErrorsPython(False)

solvers= []

if "SKIP_WORHP_TESTS" not in os.environ and ca.has_nlpsol("worhp")  and not args.ignore_memory_heavy:
  solvers.append(("worhp",{"worhp": {"TolOpti":1e-9}},{"codegen": False,"discrete":False}))
  #solvers.append(("worhp",{"TolOpti":1e-20,"TolFeas":1e-20,"UserHM": False}))
  pass

if "SKIP_FATROP_TESTS" not in os.environ and ca.has_nlpsol("fatrop"):
  flags = []
  if os.name != 'nt':
    flags = ["-Wno-strict-prototypes"]
  codegen = {"std": "c99","extralibs": ["fatrop","blasfeo"],"extra_options":flags}
  solvers.append(("fatrop",{"fatrop": {}},{"codegen":codegen,"discrete":False}))

if "SKIP_IPMC_TESTS" not in os.environ and ca.has_nlpsol("ipmc"):
  # ipmc is not installed into casadi's tree, so codegen has to be pointed
  # at its build directory; without IPMC_ROOT there is nothing to link.
  ipmc_flags = [] if os.name == 'nt' else ["-Wno-strict-prototypes"]
  ipmc_root = os.environ.get("IPMC_ROOT", "/missing")
  ipmc_libdir = os.path.join(ipmc_root, "build", "ipmc")
  ipmc_codegen = False
  if not os.path.isdir(ipmc_libdir):
    print("IPMC_ROOT not set or has no build/ipmc, skipping ipmc codegen checks")
  elif not os.path.isdir(ca.GlobalOptions.getCasadiIncludePath()):
    # The generated code pulls in ipmc's header, which pulls in <blasfeo.h>
    # from casadi's own include tree -- absent in a non-selfcontained install.
    # Same limitation applies to the fatrop/madnlp/uno codegen checks.
    print("no casadi include tree at %s, skipping ipmc codegen checks"
          % ca.GlobalOptions.getCasadiIncludePath())
  else:
    # headers live in the source tree, the library in the build tree
    ipmc_codegen = {"std": "c99","extralibs": ["ipmc","blasfeo"],
                       "extralibdirs": [ipmc_libdir],"extra_include": [ipmc_root],
                       "extra_options": ipmc_flags}
  solvers.append(("ipmc",{"ipmc": {}},{"codegen": ipmc_codegen,"discrete":False}))

if "SKIP_SLEQP_TESTS" not in os.environ and ca.has_nlpsol("sleqp"):
  solvers.append(("sleqp",{"print_time":False,"sleqp": {"linesearch": "Approx","feas_tol":1e-7,"stat_tol":1e-7,"slack_tol":1e-7, "hess_eval": "Exact"}},{"codegen": False,"discrete":False}))

if "SKIP_UNO_TESTS" not in os.environ and ca.has_nlpsol("uno"):
  uno_codegen = {"std": "c99", "extralibs": ["uno"],"extra_include": ["uno"]}
  solvers.append(("uno",{"print_time":False,"uno": {"preset": "ipopt", "primal_tolerance":1e-8, "dual_tolerance":1e-8}}, {"codegen": uno_codegen, "discrete": False}))
  solvers.append(("uno",{"print_time":False,"uno": {"preset": "funnelsqp", "primal_tolerance":1e-8, "dual_tolerance":1e-8}}, {"codegen": uno_codegen, "discrete": False}))
  solvers.append(("uno",{"print_time":False,"uno": {"preset": "filtersqp", "primal_tolerance":1e-8, "dual_tolerance":1e-8}}, {"codegen": uno_codegen, "discrete": False}))
  solvers.append(("uno",{"print_time":False,"uno": {"preset": "filtersqp", "primal_tolerance":1e-8, "dual_tolerance":1e-8, "hessian_model": "LBFGS"}}, {"codegen": uno_codegen, "discrete": False}))

if "SKIP_ALPAQA_TESTS" not in os.environ and ca.has_nlpsol("alpaqa"):
  solvers.append(("alpaqa",{"print_time":False,"alpaqa": {"alm.tolerance": 1e-10, "alm.dual_tolerance": 1e-10, "alm.penalty_update_factor": 10, "alm.max_iter": 3000, "alm.print_interval": 1, "panoc.max_iter": 500, "panoc.print_interval": 1, "lbfgs.memory": 2}},{"codegen": False,"discrete":False}))

if "SKIP_IPOPT_TESTS" not in os.environ and ca.has_nlpsol("ipopt"):
  codegen = {"extralibs": ["ipopt"], "std": "c99"}
  # casadi's generated C for ipopt is written against ipopt >= 3.14's C
  # interface, which spells the callback types ipindex / ipnumber / bool.
  # 3.11 - 3.13 spell them Index / Number / Bool, so the generated file does
  # not compile at all there (every callback signature is rejected).  Only
  # DISABLE on a positive sighting of such a header, so that a tree whose
  # header we cannot locate behaves exactly as before.
  for _h in [os.path.join(ca.GlobalOptions.getCasadiIncludePath(), "coin-or", "IpStdCInterface.h"),
             "/usr/include/coin-or/IpStdCInterface.h",
             "/usr/local/include/coin-or/IpStdCInterface.h"]:
    if os.path.isfile(_h):
      with open(_h) as _f:
        if "ipindex" not in _f.read():
          print("ipopt < 3.14 (%s does not declare 'ipindex'), "
                "skipping ipopt codegen checks" % _h)
          codegen = False
      break
  solvers.append(("ipopt",{"print_time":False,"ipopt": {"tol": 1e-10, "derivative_test":"second-order","print_level":0}},{"codegen": codegen,"discrete":False}))
  solvers.append(("ipopt",{"print_time":False,"ipopt": {"tol": 1e-10, "derivative_test":"first-order","hessian_approximation": "limited-memory","print_level":0}},{"codegen": codegen,"discrete":False}))

if "SKIP_IPOPT_TESTS" not in os.environ and ca.has_nlpsol("ipopt") and ca.has_nlpsol("sqpmethod"):
  qpsol_options = {"nlpsol": "ipopt", "nlpsol_options": {"ipopt.tol": 1e-12,"ipopt.tiny_step_tol": 1e-20, "ipopt.fixed_variable_treatment":"make_constraint","ipopt.print_level":0,"print_time":False,"print_time":False} }
  solvers.append(("sqpmethod",{"qpsol": "nlpsol","qpsol_options": qpsol_options,"print_header":False,"print_iteration":True,"print_time":False},{"codegen": False,"discrete":False}))
  solvers.append(("sqpmethod",{"qpsol": "nlpsol","qpsol_options": qpsol_options,"hessian_approximation": "limited-memory","tol_du":1e-10,"tol_pr":1e-10,"min_step_size":1e-14,"print_header":False,"print_iteration":True,"print_time":False,"max_iter":55},{"codegen": False,"discrete":False}))

if "SKIP_QRQP_TESTS" not in os.environ and ca.has_conic("qrqp") and ca.has_nlpsol("sqpmethod"):
  codegen = {"std": "c99"}
  qpsol_options = {"print_iter":False,"print_header":False,"error_on_fail" : False}
  solvers.append(("sqpmethod",{"qpsol": "qrqp","qpsol_options": qpsol_options,"print_header":False,"print_iteration":False,"print_time":False},{"codegen":codegen,"discrete":False}))
  solvers.append(("sqpmethod",{"qpsol": "qrqp","max_iter_ls":0,"qpsol_options": qpsol_options,"print_header":False,"print_iteration":False,"print_time":False},{"codegen":codegen,"discrete":False}))
  solvers.append(("sqpmethod",{"qpsol": "qrqp","convexify_strategy":"regularize","max_iter":500,"qpsol_options": qpsol_options,"print_header":False,"print_iteration":True,"print_time":False,"tol_du":1e-8,"min_step_size":1e-12},{"codegen":codegen,"discrete":False}))

if "SKIP_DAQP_TESTS" not in os.environ and ca.has_conic("daqp") and ca.has_nlpsol("sqpmethod"):
  codegen = {"std": "c99","extralibs": ["daqp"]}
  solvers.append(("sqpmethod",{"qpsol": "daqp","hessian_approximation":"limited-memory","tol_du":1e-10,"tol_pr":1e-10,"min_step_size":1e-14,"max_iter":55,"print_header":False,"print_iteration":True,"print_time":False,"qpsol_options":{"print_problem":True,"print_out":True}},{"codegen":codegen,"discrete":False}))

if "SKIP_BLOCKSQP_TESTS" not in os.environ and ca.has_nlpsol("blocksqp"):
  try:
    ca.load_linsol("ma27")
    solvers.append(("blocksqp",{},{"codegen": False,"discrete":False}))
  except:
    pass

if "SKIP_BONMIN_TESTS" not in os.environ and ca.has_nlpsol("bonmin"):
  solvers.append(("bonmin",{},{"discrete":True,"codegen":False}))

if "SKIP_KNITRO_TESTS" not in os.environ and ca.has_nlpsol("knitro"):
  solvers.append(("knitro",{"knitro":{"feastol":1e-9,"opttol":1e-9,"ftol":1e-20}},{"codegen": False,"discrete":False}))

if "SKIP_SNOPT_TESTS" not in os.environ and ca.has_nlpsol("snopt"):
  solvers.append(("snopt",{"snopt": {"Verify_level": 3,"Major_optimality_tolerance":1e-12,"Minor_feasibility_tolerance":1e-12,"Major_feasibility_tolerance":1e-12}},{"codegen": False,"discrete":False}))

libmad_dir = os.environ.get("LIBMADDIR", "/missing")
libmad_codegen = {"extralibs": ["Mad"], "std": "c99",
                "extralibdirs": [os.path.join(libmad_dir, "lib")],
                "extra_include": [os.path.join(libmad_dir, "include")]}

if "SKIP_MADNLP_TESTS" not in os.environ and ca.has_nlpsol("madnlp"):
  solvers.append(("madnlp",{"madnlp": {}},{"codegen": libmad_codegen,"discrete":False}))

print(solvers)

class NLPtests(casadiTestCase):

  @requires_nlpsol("ccopt")
  def test_ccopt(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    f = (x-2)**2 + (y-2)**2
    nlp = {"x": ca.vertcat(x, y), "f": f}

    solver_opts = {
      "cc_pairs": [[0, 1]],
      "madnlp": {"bound_relax_factor": 0.0},
    }
    solver = ca.nlpsol("solver", "ccopt", nlp, solver_opts)

    solver_in = {
      "x0": [2.5, 3.0],
      "lbx": [0.0, 0.0],
    }
    solver_out = solver(**solver_in)

    self.assertTrue(solver.stats()["success"])
    self.assertAlmostEqual(float(solver_out["f"]), 4.0, 6)
    x_sol = float(solver_out["x"][0])
    y_sol = float(solver_out["x"][1])
    self.assertAlmostEqual(min(x_sol, y_sol), 0.0, 6)
    self.assertAlmostEqual(max(x_sol, y_sol), 2.0, 6)

    aux_codegen = {"extralibs": ["Mad"], "std": "c99",
                   "extralibdirs": [os.path.join(os.environ.get("LIBMADDIR", "/missing"), "lib")],
                   "extra_include": [os.path.join(os.environ.get("LIBMADDIR", "/missing"), "include")]}
    self.check_codegen(solver, solver_in, **aux_codegen)
    self.check_serialize(solver, solver_in)

  @requires_nlpsol("alpaqa")
  def test_alpaqa(self):
  
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}

    solver = ca.nlpsol("mysolver", "alpaqa", nlp)
    solver_in = {}
    solver_in["lbx"]=[-10]*2
    solver_in["ubx"]=[10]*2
    solver_in["lbg"]=[-10]
    solver_in["ubg"]=[10]
    solver_out = solver(**solver_in)

  @memory_heavy()
  def test_nonregular_point(self):
    x=ca.SX.sym("x")

    nlp={'x':x,'f':(x+1)**2, 'g': ca.sqrt(x)}

    for Solver, solver_options, aux_options in solvers:
      # ipmc loops forever instead of bailing out on the NaN
      if Solver in ["madnlp", "fatrop", "ipmc"]: continue

      print("test_nonregular_point",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-1000]
      solver_in["ubx"]=[1000]
      solver_in["lbg"]=[-1000]
      solver_in["ubg"]=[1000]
      solver_in["x0"] = 1e-4
      try:
        print(solver(**solver_in))
      except:
        pass
      if Solver not in ["ipopt","snopt","blocksqp","bonmin","knitro","sleqp","alpaqa","uno"]:
        self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_NAN")
      self.assertFalse(solver.stats()["success"])

    nlp={'x':x,'f':x**2, 'g': ca.sqrt(x)}
    for Solver, solver_options, aux_options in solvers:
      # ipmc loops forever instead of bailing out on the NaN
      if Solver in ["madnlp", "fatrop", "ipmc"]: continue
      print("test_nonregular_point",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-2]
      solver_in["ubg"]=[2]
      solver_in["x0"] = -2
      try:
        print(solver(**solver_in))
        solver_out = solver(**solver_in)
      except:
        pass
      if Solver not in ["ipopt","snopt","bonmin","knitro","sleqp","alpaqa"]:
        self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_NAN")
      self.assertFalse(solver.stats()["success"])

  def test_iteration_interrupt(self):
  
   #add test for keyboard interrupt in fun_eval, not just iteration_callback
   
   class mycallback(ca.Callback):
      def __init__(self, name):
        ca.Callback.__init__(self)
        opts = {"enable_fd":True}
        self.construct(name, opts)

      def eval(self,argin):
        raise KeyboardInterrupt()

   interrupt = mycallback("interrupt")  # pyright: ignore[reportAssignmentType]
   x = ca.MX.sym("x")
   nlp = {"x":x,"f":interrupt(x),"g":x}
   for Solver, solver_options, aux_options in solvers:
     if Solver in ["madnlp"]: continue
     solver_options = dict(solver_options)
     solver_options["error_on_fail"] = True
     solver = ca.nlpsol("solver",Solver,nlp,solver_options)

     with self.assertInAnyOutput("KeyboardInterrupt"):
       solver(lbg=-5,ubg=5)

     with self.assertRaises(Exception):
         solver(lbg=-5,ubg=5)

   for Solver, solver_options, aux_options in solvers:
      
      #if Solver not in ["ipopt","sqpmethod"]: continue
      if Solver in ["worhp","blocksqp","knitro","bonmin","snopt","alpaqa","fatrop","madnlp"]: continue
      print("test_iteration_interrupt",Solver,solver_options)

      opti = ca.Opti()

      x = opti.variable()

      def interrupt(i):
        raise KeyboardInterrupt()

      opti.minimize((x-1)**4)
      opti.subject_to(-5<=(x<=5))


      opts = dict(solver_options)
      if Solver=="bonmin":
        opts["discrete"] = [1]

      opti.solver(Solver, opts)
      sol = opti.solve()
      opti.callback(interrupt)
      with self.assertRaises(Exception):
        sol = opti.solve()

      def interrupt(i):
        raise Exception()
      opti.callback(interrupt)
      opts["iteration_callback_ignore_errors"] = False
      opti.solver(Solver, opts)
      with self.assertRaises(Exception):
        sol = opti.solve()

      opts["iteration_callback_ignore_errors"] = True
      opti.solver(Solver, opts)
      sol = opti.solve()

      def interrupt(i):
        raise KeyboardInterrupt()
      opti.callback(interrupt)
      with self.assertRaises(Exception):
        sol = opti.solve()


  def test_discrete(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1.4-x)**2+100*(y-x**2)**2}

    for Solver, solver_options, aux_options in solvers:
      if not aux_options["discrete"]: continue
      print("test_discrete",Solver,solver_options)
      solver_options = dict(solver_options)
      solver_options["discrete"] = [1,0]
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0.16,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))

      self.check_serialize(solver,solver_in)

    for Solver, solver_options, aux_options in solvers:
      if not aux_options["discrete"]: continue
      print("test_discrete",Solver,solver_options)
      solver_options = dict(solver_options)
      solver_options["discrete"] = [0,1]
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],2e-4,3,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],ca.sqrt(2),3,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],2,3,str(Solver))

      self.check_serialize(solver,solver_in)
  def test_nan(self):
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':-x,'g':x}
    for Solver, solver_options, aux_options in solvers:
      
      print("test_nan",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)

      for x in ["x","g"]:
        lb = "lb"+x
        ub = "ub"+x
        for data in [{lb:3,ub:-3},
                     {lb:np.inf,ub:np.inf},
                     {lb:-np.inf,ub:-np.inf},
                     {lb:np.nan},
                     {ub:np.nan},
                     ]:
          print(data)
          with self.assertInException("Ill-posed"):
            solver(**data)

  # ---------------------------------------------------------------------
  # Soft constraints (the 'S' option / 's' / 'f_s')
  #
  # The unconstrained optimum sits at (3, 2). g[0] and the upper bound on
  # x[0] push down on it while g[2] pushes up, so both the s_l and the s_u
  # half of the slack vector end up active. g[1] and the bounds on x[1] are
  # left hard, which exercises the hard/soft reordering of the constraint
  # vector too.
  #
  # See docs/examples/python/nlpsol_slacks.py for a worked example, and
  # nlpsol_slacks_manual.py for the augmentation done by hand.
  # ---------------------------------------------------------------------
  def slack_problem(self):
    x = ca.SX.sym("x", 2)
    f = (x[0]-3)**2 + (x[1]-2)**2
    g = ca.vertcat(x[0]+x[1], x[0]-x[1], x[0]+2*x[1])
    bounds = {"x0": [0, 0], "lbx": [-10, -10], "ubx": [0.5, 10],
              "lbg": [-inf, -inf, 4.0], "ubg": [1.0, 0.5, inf]}
    return x, f, g, bounds

  def slack_norms(self):
    w = 1.0
    # g[0] soft, g[1] hard, g[2] soft, x[0] soft, x[1] hard
    S_sep = ca.sparsify(ca.DM([[1, 0, 0],
                               [0, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1],
                               [0, 0, 0]])).sparsity()
    # a single slack shared by all softened rows
    S_shared = ca.sparsify(ca.DM([[1], [0], [1], [1], [0]])).sparsity()
    return [("L1",   S_sep,    lambda s: w*ca.sum1(s)),
            ("L2",   S_sep,    lambda s: w*ca.dot(s, s)),
            ("Linf", S_shared, lambda s: w*ca.sum1(s))]

  def slack_reference(self, Solver, solver_options, S, penalty,
                      par=None, pval=None):
    """Hand-augmented reference: what nlpsol's slack layer does for you.

    Deliberately written differently from the internal de-sugaring: every row
    is relaxed here (the empty rows of S simply reproduce the original bound),
    so agreement is not an artefact of a shared layout.  And -- where ipopt is
    available -- deliberately solved by a DIFFERENT solver, so agreement is
    not an artefact of a shared solver core either.  'Solver'/'solver_options'
    are then used only as the fallback for a tree without ipopt.

    'par' is an optional parameter the penalty depends on; 'pval' its value
    for this solve (see test_slacks_parametric_penalty)."""
    x, f, g, bounds = self.slack_problem()
    nx, ng = x.numel(), g.numel()
    ns = S.size2()

    Sd = ca.DM.ones(S)  # S is structural; its entries count as 1
    S_g, S_x = Sd[:ng, :], Sd[ng:, :]

    s = ca.SX.sym("s", 2*ns)
    s_l, s_u = s[:ns], s[ns:]

    G = ca.vertcat(g + ca.mtimes(S_g, s_l),   # >= lbg
                   g - ca.mtimes(S_g, s_u),   # <= ubg
                   x + ca.mtimes(S_x, s_l),   # >= lbx
                   x - ca.mtimes(S_x, s_u))   # <= ubx
    lbG = ca.vertcat(ca.DM(bounds["lbg"]), -inf*ca.DM.ones(ng),
                     ca.DM(bounds["lbx"]), -inf*ca.DM.ones(nx))
    ubG = ca.vertcat(inf*ca.DM.ones(ng), ca.DM(bounds["ubg"]),
                     inf*ca.DM.ones(nx), ca.DM(bounds["ubx"]))

    rnlp = {"x": ca.vertcat(x, s), "f": f + penalty(s), "g": G}
    extra = {}
    if par is not None:
      rnlp["p"] = par
      extra["p"] = pval
    # Solved by IPOPT, not by the solver under test: see the note below
    # slack_reference() for why, and for what happens without ipopt.
    if self.slack_ref_solver is None:
      RefSolver, ref_options = Solver, solver_options
    else:
      RefSolver = self.slack_ref_solver
      ref_options = dict(self.slack_ref_opts,
                         ipopt=dict(self.slack_ref_opts["ipopt"]))
    solver = ca.nlpsol("reference", RefSolver, rnlp, ref_options)
    r = solver(x0=ca.vertcat(ca.DM(bounds["x0"]), ca.DM.zeros(2*ns)),
               lbx=ca.vertcat(-inf*ca.DM.ones(nx), ca.DM.zeros(2*ns)),
               ubx=inf*ca.DM.ones(nx+2*ns), lbg=lbG, ubg=ubG, **extra)

    # Map the canonical solution back onto the user-facing quantities
    z, lamZ, lamG = r["x"], r["lam_x"], r["lam_g"]
    return {"f": r["f"],
            "x": z[:nx],
            "s": z[nx:],
            "g": ca.Function("g", [x], [g])(z[:nx]),
            "lam_s": lamZ[nx:],
            # at most one of the two one-sided rows of a relaxed row is active
            "lam_g": lamG[:ng] + lamG[ng:2*ng],
            "lam_x": lamZ[:nx] + lamG[2*ng:2*ng+nx] + lamG[2*ng+nx:]}

  # A solver that de-sugars slacks internally (sqpmethod and friends) solves
  # *literally the same NLP* on both sides of the comparison below, so it
  # agrees with slack_reference() to solver precision at any tolerance.  A
  # solver with a NATIVE slack path (ipmc) solves two genuinely different
  # formulations of the same problem, which only have to agree at the optimum
  # -- and the optimum of slack_problem() is a degenerate vertex: at the L1
  # solution x = [2, 1.5] the hard row g[1] = x[0]-x[1] <= 0.5 is active with
  # a multiplier of exactly zero, so strict complementarity fails.  An
  # interior-point method approaches such a point only at O(sqrt(mu)), along a
  # central path that depends on the formulation, and at ipmc's default tol
  # the two formulations land ~3e-6 apart in x -- both being about that far
  # from the true optimum.  Against an independent high-accuracy sqpmethod
  # solve the native answer is in fact *closer* than the reference
  # augmentation (2.8e-5 vs 3.1e-5 at tol 1e-8; both ~2e-7 at tol 1e-12), so
  # the disagreement is a convergence artefact of the shared degenerate
  # vertex, not a slack-mapping error.
  #
  # The fix is therefore to solve *both* sides more accurately (tol=1e-13),
  # not to weaken the comparison.
  #
  # The reference itself is now IPOPT rather than the solver under test
  # (slack_ref_solver below).  With the solver under test on both sides its
  # own convergence error cancelled exactly, so the comparison could only
  # ever see a slack-MAPPING mistake, never a wrong answer that both sides
  # shared; against ipopt it cannot cancel.  The price is that each solver's
  # own O(sqrt(mu)) approach to that degenerate vertex is now visible, which
  # is what sets the tolerances below.  Where ipopt is absent the reference
  # falls back to the solver under test and the test keeps its original
  # (weaker, but not vacuous) meaning.
  #
  # Reference options: tol 1e-14 with mu_strategy=adaptive, and
  # bound_relax_factor=0.  Both matter and both were measured.  At tol 1e-12
  # the reference is itself 4.7e-7 away from the exact L1 vertex (checked
  # against sqpmethod/qrqp, an active-set method with no barrier at all);
  # tol 1e-14 + adaptive brings that to 3.8e-8.  bound_relax_factor=0 stops
  # the reference from widening every bound by 1e-8 the way ipopt and
  # ipmc both do by default, which is worth 6e-8 on f here.
  #
  # ipopt is ALSO in `solvers`, and at its listed tol of 1e-10 it lands
  # 3.3e-6 from the vertex -- ten times the worst of any other solver.  It is
  # therefore tightened for these tests exactly like ipmc is; that is a
  # tightening of the solve, not a weakening of the check.
  #
  # Tolerances, measured over every solver in `solvers` x every norm, against
  # the ipopt reference: f 6.9e-8, x 4.3e-7, g 8.7e-7, lam_g 8.7e-7,
  # lam_x 8.0e-9 (worst three of those at sqpmethod/L1).  Hence digits=6 on
  # f and digits=5 on the rest -- x, s and g drop a digit from the era when
  # the reference was the solver under test, and that drop is the visible
  # part of the cancellation described above, not slack-layer error.
  #
  # The exception is again the L2 case's should-be-zero slacks: an L2 penalty
  # has zero gradient at s = 0, so a slack that ought to be zero stalls on the
  # interior-point barrier floor sqrt(mu/2Z) of whichever formulation it lives
  # in, at a slightly different mu.  For the native path that is s 2.9e-6 and
  # lam_s 5.7e-6 (against 6.4e-7 / 1.3e-6 for the de-sugaring solvers), so
  # digits=4 on 's'/'lam_s' -- for the native path only, for the L2 norm only.
  slack_native = ["ipmc"]           # solvers with a NATIVE slack path
  slack_tight_tol = {"ipmc": 1e-13, "ipopt": 1e-13}
  # Which solvers the slack tests run.  An allow-list rather than a skip list:
  # these tests validate the slack LAYER, and the problem they use is
  # deliberately degenerate (at the optimum the hard row g[1] is active with a
  # multiplier of exactly zero, so strict complementarity fails).  Every solver
  # approaches such a point on its own central path and at its own default
  # tolerance, so "all solvers minus the ones that happened to fail" grew a new
  # entry each time one was excluded -- fatrop stalls at max_iter, sleqp refuses
  # the expanded problem outright, uno lands 3e-5 away in x and does not report
  # success on the parametric penalty.  None of that is about the slack layer.
  # ipmc drives the NATIVE path, ipopt is the reference the others are compared
  # against, and sqpmethod is an independent de-sugaring solver: between them
  # both paths through the layer are covered.
  slack_solvers = ["ipopt", "ipmc", "sqpmethod"]
  slack_ref_solver = "ipopt" if ca.has_nlpsol("ipopt") else None
  slack_ref_opts = {"print_time": False,
                    "ipopt": {"print_level": 0, "sb": "yes", "tol": 1e-14,
                              "mu_strategy": "adaptive",
                              "constr_viol_tol": 1e-14, "compl_inf_tol": 1e-14,
                              "dual_inf_tol": 1e-9, "acceptable_tol": 1e-14,
                              "bound_relax_factor": 0, "max_iter": 5000}}

  def test_slacks(self):
    x, f, g, bounds = self.slack_problem()

    for Solver, base_options, aux_options in solvers:
      if Solver not in self.slack_solvers: continue
      # sqpmethod with an L-BFGS Hessian never reports success on this
      # problem, whatever it is asked to solve: the plain hard NLP with no
      # slacks at all stops after 3 iterations with success=False, and so do
      # the native slack path (55 iterations, its max_iter) and the hand
      # augmentation.  Its answers are right to ~5e-7 in every case, but
      # there is no converged solve to assert against, so it is skipped here
      # the same way fatrop is.  Nothing to do with the slack layer; note
      # that this configuration only exists at all when ipopt is built (it
      # uses ipopt as its qpsol).  test_slacks_ubs still exercises it.
      if Solver == "sqpmethod" and \
         base_options.get("hessian_approximation") == "limited-memory":
        continue
      solver_options = dict(base_options)
      if Solver in self.slack_tight_tol:
        solver_options[Solver] = dict(solver_options.get(Solver, {}),
                                      tol=self.slack_tight_tol[Solver])
      for name, S, penalty in self.slack_norms():
        print("test_slacks", Solver, name, solver_options)
        # See the note above: only the native-slack path on the L2 penalty
        # loosens, and only on the slacks that should be zero.
        digits = {"f": 6, "x": 5, "s": 5, "g": 5,
                  "lam_x": 5, "lam_g": 5, "lam_s": 5}
        if Solver in self.slack_native and name == "L2":
          digits["s"] = 4
          digits["lam_s"] = 4
        ns = S.size2()
        s = ca.SX.sym("s", 2*ns)
        nlp = {"x": x, "f": f, "g": g, "s": s, "f_s": penalty(s)}
        opts = dict(solver_options)
        opts["S"] = S

        solver = ca.nlpsol("mysolver", Solver, nlp, opts)
        self.assertEqual(solver.n_in(), ca.nlpsol_n_in())
        self.assertEqual(solver.size1_in("ubs"), 2*ns)
        self.assertEqual(solver.size1_out("s"), 2*ns)
        # 'x'/'g' keep the sizes the user wrote, not the augmented ones
        self.assertEqual(solver.size1_out("x"), x.numel())
        self.assertEqual(solver.size1_out("g"), g.numel())

        r = solver(**bounds)
        self.assertTrue(solver.stats()["success"])

        ref = self.slack_reference(Solver, solver_options, S, penalty)
        for k in ["f", "x", "s", "g", "lam_x", "lam_g", "lam_s"]:
          self.checkarray(r[k], ref[k], name+":"+k, digits=digits[k])

        # The solution is feasible for the relaxed bounds, hard rows included
        Sd = ca.DM.ones(S)
        z = ca.vertcat(r["g"], r["x"])
        lb = ca.vertcat(ca.DM(bounds["lbg"]), ca.DM(bounds["lbx"]))
        ub = ca.vertcat(ca.DM(bounds["ubg"]), ca.DM(bounds["ubx"]))
        self.assertTrue(float(ca.mmax(z - ub - ca.mtimes(Sd, r["s"][ns:]))) < 1e-6)
        self.assertTrue(float(ca.mmax(lb - ca.mtimes(Sd, r["s"][:ns]) - z)) < 1e-6)

        if aux_options["codegen"]:
          self.check_codegen(solver, bounds, **aux_options["codegen"])
        self.check_serialize(solver, bounds)

  def test_slacks_ubs(self):
    """ubs=0 forbids any relaxation, so the problem must reduce to the hard NLP."""
    x, f, g, bounds = self.slack_problem()

    for Solver, solver_options, aux_options in solvers:
      if Solver not in self.slack_solvers: continue
      for name, S, penalty in self.slack_norms():
        print("test_slacks_ubs", Solver, name, solver_options)
        ns = S.size2()
        s = ca.SX.sym("s", 2*ns)
        opts = dict(solver_options)
        opts["S"] = S

        soft = ca.nlpsol("mysolver", Solver, {"x": x, "f": f, "g": g,
                                              "s": s, "f_s": penalty(s)}, opts)
        hard = ca.nlpsol("hard", Solver, {"x": x, "f": f, "g": g}, solver_options)

        r = soft(ubs=ca.DM.zeros(2*ns), **bounds)
        rh = hard(**bounds)
        self.checkarray(r["s"], ca.DM.zeros(2*ns), name+":s", digits=6)
        self.checkarray(r["x"], rh["x"], name+":x", digits=6)
        self.checkarray(r["f"], rh["f"], name+":f", digits=6)
        self.checkarray(r["g"], rh["g"], name+":g", digits=6)

  def test_slacks_parametric_penalty(self):
    """f_s depends on the parameter p, and weights the two halves of s
    differently.

    Two gaps in one problem.

    (a) A p-dependent penalty is the one case a native slack path can
    silently get wrong, by reading Z and z once and then solving every
    later call with the weights p happened to hold then.  Nothing else in
    the suite exercises that, so an implementation that evaluated the
    penalty once at construction, or constant-folded it into the generated
    code, would pass everything else.  ONE solver object is therefore
    asked for three solves, and the third repeats the first: a stale cache
    cannot reproduce it.  (ipmc is excluded below: it does bake z and Z in
    at construction, and refuses a p-dependent f_s outright rather than
    solving it wrongly -- see ocp.py's
    test_ipmc_slacks_parametric_penalty_refused.)

    (b) Every other slack test weights s_l and s_u identically, which
    makes a zl/zu (or S_l/S_u) swap invisible.  Here they cost 0.45 and
    2.5, and the case is run both ways round.  At (0.45, 2.5) both halves
    are active at once -- s_l = [0, 0.875, 0] and s_u = [0.925, 0, 0.225]
    -- and swapping the weights moves the optimum from f = 9.08 to
    f = 2.71, so the two are in no way interchangeable."""
    x, f, g, bounds = self.slack_problem()
    S = self.slack_norms()[0][1]           # the L1 sparsity: g[0], g[2], x[0]
    ns = S.size2()
    s = ca.SX.sym("s", 2*ns)
    par = ca.SX.sym("w_slack", 2)
    penalty = lambda s: par[0]*ca.sum1(s[:ns]) + par[1]*ca.sum1(s[ns:])
    weights = [ca.DM([0.45, 2.5]), ca.DM([2.5, 0.45])]

    for Solver, base_options, aux_options in solvers:
      if Solver not in self.slack_solvers: continue
      # see test_slacks: sqpmethod/L-BFGS never reports success here
      if Solver == "sqpmethod" and \
         base_options.get("hessian_approximation") == "limited-memory":
        continue
      # ipmc hands ipmc plain numbers for z and Z, evaluated once at
      # construction, so it refuses a p-dependent penalty instead of
      # solving it with stale weights.  The refusal, and the
      # expand_slacks route that does track p, are tested in ocp.py's
      # test_ipmc_slacks_parametric_penalty_refused.
      if Solver == "ipmc": continue
      solver_options = dict(base_options)
      if Solver in self.slack_tight_tol:
        solver_options[Solver] = dict(solver_options.get(Solver, {}),
                                      tol=self.slack_tight_tol[Solver])
      print("test_slacks_parametric_penalty", Solver, solver_options)
      opts = dict(solver_options)
      opts["S"] = S
      solver = ca.nlpsol("mysolver", Solver,
                         {"x": x, "f": f, "g": g, "p": par,
                          "s": s, "f_s": penalty(s)}, opts)

      out = []
      for pv in [weights[0], weights[1], weights[0]]:
        r = solver(p=pv, **bounds)
        self.assertTrue(solver.stats()["success"])
        out.append(r)
        ref = self.slack_reference(Solver, solver_options, S, penalty,
                                   par=par, pval=pv)
        for k in ["f", "x", "s", "g"]:
          self.checkarray(r[k], ref[k],
                          "par%s:%s" % (str(pv), k), digits=6)

      # Non-vacuity. At the first weighting BOTH halves of s are active,
      # so the test can see which weight was applied where; the swapped
      # weighting is a materially different answer, so a solve that
      # ignored p could not match both.
      sl = [float(out[0]["s"][i]) for i in range(ns)]
      su = [float(out[0]["s"][ns+i]) for i in range(ns)]
      n_l = sum(1 for v in sl if v > 1e-4)
      n_u = sum(1 for v in su if v > 1e-4)
      print("test_slacks_parametric_penalty", Solver, "active", (n_l, n_u),
            "f", float(out[0]["f"]), float(out[1]["f"]))
      self.assertTrue(n_l >= 1 and n_u >= 2)
      self.assertTrue(abs(float(out[0]["f"]-out[1]["f"])) > 1.0)
      self.assertTrue(float(ca.norm_inf(out[0]["x"]-out[1]["x"])) > 0.5)
      # Same p, same solver object, same answer -- nothing from the middle
      # solve may leak into the third.
      for k in ["f", "x", "s", "g"]:
        self.checkarray(out[0][k], out[2][k], "par_repeat:"+k, digits=10)

      if aux_options["codegen"]:
        # the generated code has to read the penalty's z from p as well
        self.check_codegen(solver, dict(bounds, p=weights[1]),
                           **aux_options["codegen"])
      self.check_serialize(solver, dict(bounds, p=weights[0]))

  @requires_nlpsol("ipmc")
  def test_slacks_s0(self):
    """s0 must reach the solver, and must not change the answer.

    The softened row here has to be relaxed a long way before anything is
    feasible: x0*x1 = 1 with x0,x1 >= 0.05 forces x0+x1 >= 2, while ubg[1] is
    -1, so s_u has to reach 3.04 and the softened row is active at the
    solution.  Where the slack starts therefore decides how the solve goes,
    and restoration is involved on the way.

    The native slack path used to discard s0 outright -- it produced
    bit-identical output for every guess a user could write down -- and once
    it stopped doing that, a cached solver object still only honoured the s0
    of its first solve, because ipmc reads the initial slack once at
    create.  Both are asserted here: the same solver object is re-solved with
    three different s0.
    """
    if "SKIP_IPMC_TESTS" in os.environ: return
    aux_options = [a for (S, _b, a) in solvers if S == "ipmc"][0]

    x = ca.SX.sym("x", 2)
    s = ca.SX.sym("s", 2)
    nlp = {"x": x, "f": (x[0]-4)**2 + x[1]**2,
           "g": ca.vertcat(x[0]*x[1], x[0]+x[1]),
           "s": s, "f_s": 20*ca.sum1(s)}
    # softens g[1] only, on its upper side and its lower side
    S = ca.Sparsity.triplet(4, 1, [1], [0])
    bounds = dict(lbg=[1.0, -ca.inf], ubg=[1.0, -1.0],
                  lbx=[0.05, 0.05], ubx=[10, 10], x0=[0.0, 0.0])
    # [0,0] is the default s0, and is the guess the solve is hardest from: the
    # relaxation has to be found from scratch, through restoration
    guesses = [[0, 0.0], [0, 2.0], [0, 3.5], [0, 6.0]]

    ref = None
    for tag, extra in [("native", {}), ("expand", {"expand_slacks": True})]:
      opts = dict({"S": S, "print_time": False,
                   "ipmc": {"tol": 1e-10, "max_iter": 300}}, **extra)
      solver = ca.nlpsol("mysolver", "ipmc", nlp, opts)
      counts = []
      for s0 in guesses:
        r = solver(s0=s0, **bounds)
        self.assertTrue(solver.stats()["success"],
                        "%s s0=%s did not converge" % (tag, s0))
        counts.append(solver.stats()["iter_count"])
        if ref is None:
          ref = dict((k, r[k]) for k in ["f", "x", "s"])
          # the softened row is genuinely active, or this proves nothing
          self.assertTrue(float(r["s"][1]) > 2.5)
        # native and expansion, and every s0, agree on the answer
        for k in ["f", "x", "s"]:
          self.checkarray(r[k], ref[k], "%s s0=%s :%s" % (tag, s0, k), digits=5)
      # s0 reaches the solver: a different starting slack takes a different
      # route.  Asserted on ONE solver object, so a cached solver that only
      # honours the first s0 fails here.
      print("test_slacks_s0", tag, counts)
      self.assertTrue(len(set(counts)) > 1,
                      "%s ignores s0, iter_count %s" % (tag, counts))

      if aux_options["codegen"]:
        self.check_codegen(solver, dict(bounds, s0=guesses[-1]),
                           **aux_options["codegen"])
      self.check_serialize(solver, dict(bounds, s0=guesses[-1]))

  @requires_nlpsol("ipmc")
  def test_infeasible_not_success(self):
    """An infeasible problem must not come back as a solution.

    ipmc reaches STATE_CONVERGED from the restoration phase as well as from
    the normal one -- restoration converging to a point that is still
    infeasible for the ORIGINAL problem is ipopt's Restoration_Failed -- and
    the two were reported identically.  The result was success=True on a point
    violating a hard equality by 0.46.

    There are no slacks anywhere in this problem; the softened twin below only
    checks that a relaxation which cannot reach feasibility is reported the
    same way.
    """
    if "SKIP_IPMC_TESTS" in os.environ: return
    x = ca.SX.sym("x", 2)
    f = (x[0]-4)**2 + x[1]**2
    g = ca.vertcat(x[0]*x[1], x[0]+x[1])
    # x0*x1 = 1 with x0,x1 >= 0.05 forces x0+x1 >= 2, so x0+x1 <= 1.5 cannot hold
    bounds = dict(lbg=[1.0, -ca.inf], ubg=[1.0, 1.5],
                  lbx=[0.05, 0.05], ubx=[10, 10], x0=[0.0, 0.0])
    opts = {"print_time": False, "ipmc": {"tol": 1e-10, "max_iter": 300}}

    solver = ca.nlpsol("mysolver", "ipmc", {"x": x, "f": f, "g": g}, opts)
    r = solver(**bounds)
    self.assertFalse(solver.stats()["success"])
    # and the point really is infeasible, i.e. the assertion above is not
    # passing for some unrelated reason
    self.assertTrue(abs(float(r["g"][0]) - 1.0) > 1e-3)

    # the same problem with a slack capped far below the 0.5 it would need
    # (both sides get the same cap: ipmc refuses a column pinned on one
    # side only, see test_slacks_errors)
    s = ca.SX.sym("s", 2)
    S = ca.Sparsity.triplet(4, 1, [1], [0])
    soft = ca.nlpsol("mysolver", "ipmc",
                     {"x": x, "f": f, "g": g, "s": s, "f_s": 20*ca.sum1(s)},
                     dict(opts, S=S))
    r = soft(ubs=[0.01, 0.01], **bounds)
    self.assertFalse(soft.stats()["success"])
    self.assertTrue(abs(float(r["g"][0]) - 1.0) > 1e-3)

  @requires_nlpsol("sqpmethod")
  def test_slacks_errors(self):
    x, f, g, bounds = self.slack_problem()
    s = ca.SX.sym("s", 6)
    S = self.slack_norms()[0][1]

    with self.assertInException("must not depend on 's'"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f+s[0], "g": g,
                                          "s": s, "f_s": ca.sum1(s)}, {"S": S})
    with self.assertInException("must not depend on 'x'"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g,
                                          "s": s, "f_s": ca.sum1(s)*x[0]}, {"S": S})
    with self.assertInException("even number of rows"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g,
                                          "s": ca.SX.sym("s", 5), "f_s": 0}, {"S": S})
    with self.assertInException("5-by-3"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g,
                                          "s": s, "f_s": ca.sum1(s)},
                {"S": ca.Sparsity.dense(4, 3)})
    with self.assertInException("'s' is required"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g}, {"S": S})
    with self.assertInException("'S' is required"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g,
                                          "s": s, "f_s": ca.sum1(s)}, {})
    with self.assertInException("cannot (yet) be combined with slacks"):
      ca.nlpsol("mysolver", "sqpmethod", {"x": x, "f": f, "g": g,
                                          "s": s, "f_s": ca.sum1(s)},
                {"S": S, "detect_simple_bounds": True})

  def test_wrongdims(self):
    x=ca.SX.sym("x",2)
    nlp={'x':x, 'f':-x[0],'g':ca.diag(x)}

    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      with self.assertInException("dense vector"):
        solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
    nlp={'x':x, 'f':-x[0],'g':x @ x.T}

    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      with self.assertInException("dense vector"):
        solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'f':ca.SX(1,1),'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      with self.assertInException("dense"):
        solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'f':ca.SX.zeros(0,0),'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'g':x}
    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)

    x = ca.vec(ca.diag(ca.SX.sym("x",2)))
    nlp={'x':x, 'f':x.T @ x,'g':x[0]}
    for Solver, solver_options, aux_options in solvers:
      print("test_wrongdims",Solver,solver_options)
      with self.assertInException("dense vector"):
        solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)


  def test_initialcond(self):
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':-ca.cos(x),'g':x}

    for Solver, solver_options, aux_options in solvers:
      if Solver=="fatrop": continue
      print("test_initialcond",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[6*pi+0.01]
      solver_in["lbx"]=-inf
      solver_in["ubx"]=inf
      solver_in["lbg"]=-100
      solver_in["ubg"]=100
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["x"][0],6*pi,6,str(Solver))

  def test_boundsviol(self):
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_boundsviol",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[-20]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      with self.assertRaises(Exception):
        solver_out = solver(**solver_in)

    for Solver, solver_options, aux_options in solvers:
      print("test_boundsviol",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[-20]
      with self.assertRaises(Exception):
        solver_out = solver(**solver_in)

  def test_IPOPT(self):
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPT",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["g"][0],1,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))
      print("unified_return_status",solver.stats()["unified_return_status"])      
      self.assertEqual(solver.stats()["unified_return_status"],"SOLVER_RET_SUCCESS")


      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,debug_mode=True,**aux_options["codegen"])

  def test_IPOPT_par(self):
    x=ca.SX.sym("x")
    p=ca.SX.sym("p")
    nlp={'x':x, 'p':p, 'f':(x-p)**2, 'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPT_par",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_in["p"]=1
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTinf",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-inf]
      solver_in["ubx"]=[inf]
      solver_in["lbg"]=[-inf]
      solver_in["ubg"]=[inf]

      if Solver in ["worhp"]:
        with self.assertRaises(Exception):
          solver_out = solver(**solver_in)
        return




      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver) + str(solver_out["x"][0]-1))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      if Solver in ["worhp"]:
        with self.assertRaises(Exception):
          solver_out = solver(**solver_in)
        return

      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver) + str(solver_out["x"][0]-1))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])



  def test_IPOPTrb(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2}

    for Solver, solver_options, aux_options in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      print("test_IPOPTrb",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTrb2(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      print("test_IPOPTrb2",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      if "knitro" not in str(Solver): self.assertTrue(solver.stats()["success"])

      digits = 6

      self.assertAlmostEqual(solver_out["f"][0],0,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,5,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
        
  def test_IPOPTrbf(self):
    self.message("rosenbrock fixed, limited-memory hessian approx")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTrbf",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      if "qrqp" in str(solver_options): solver_in["x0"]=[0.6,1]
      solver_in["lbx"]=[-10,1]
      solver_in["ubx"]=[10,1]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]

      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,7,str(Solver))
      if "stabilizedsqp" not in str(Solver):
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,5,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  @requires_nlpsol("knitro")
  def test_knitro_options_file(self):
    x = ca.MX.sym("x")
    
    nlp = {"x": x, "f":x**2}
    with open("my_knitro.opt","w") as f:
        f.write("maxit            42\n")
    
    solver = ca.nlpsol("solver","knitro",nlp,{"options_file":"my_knitro.opt"})
    with self.assertOutput(["maxit:","42"],[]):
        solver()
 
  @requires_nlpsol("ipopt")
  def test_ipopt_linear_solvers(self):

    solvers = ["mumps"]
    if 'spral' in ca.CasadiMeta.feature_list():
        solvers.append("spral")
    # mac's mockup-HSL ma27 is non-functional for ipopt (works on Linux)
    if ca.has_linsol("ma27") and sys.platform != "darwin":
        solvers.append("ma27")

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]
    
    for solver in solvers:
      print("test_ipopt_solver",solver)
      solver = ca.nlpsol("mysolver", "ipopt", nlp, {"ipopt.linear_solver": solver})
      solver_in = {}
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_out = solver(**solver_in)
      solver_out = solver_out
      digits = 7
      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(solver))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(solver))
      print(solver.stats())
    
  @requires_nlpsol("uno")
  def test_uno_linear_solvers(self):

    solvers = ["MUMPS"]
    if 'spral' in ca.CasadiMeta.feature_list():
        solvers.append("SSIDS")
    if ca.has_linsol("ma27") and sys.platform != "darwin":
        solvers.append("MA27")

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    for solver in solvers:
      print("test_uno_linear_solvers",solver)
      mysolver = ca.nlpsol("mysolver", "uno", nlp, {"print_time":False,
        "uno": {"preset": "ipopt", "linear_solver": solver, "primal_tolerance":1e-8, "dual_tolerance":1e-8}})
      solver_in = {}
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_out = mysolver(**solver_in)
      digits = 6
      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(solver))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(solver))
      print(mysolver.stats())

  @requires_nlpsol("ipopt")
  def test_ipopt_stats_per_call(self):
    x = ca.MX.sym("x")

    nlp = {"x":x,"f":x**2}

    solver = ca.nlpsol("solver","ipopt",nlp)

    solver()
    s1 = solver.stats()
    solver()

    s2 = solver.stats()
    self.assertEqual(s1["n_call_nlp_f"],s2["n_call_nlp_f"])

  def test_warmstart(self):

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    
    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    for Solver, solver_options, aux_options in solvers:
      print("test_warmstart",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_out = solver(**solver_in)
      oldsolver_out = solver_out

      digits = 5

      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5 if Solver=="snopt" else 8,str(Solver))
      if Solver not in ["bonmin","sleqp"]: self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,digits=codegen_check_digits,**aux_options["codegen"])

      self.message(":warmstart")
      if "ipopt" in str(Solver):
        oldsolver=solver
        options = dict(solver_options)
        options["ipopt.warm_start_init_point"]="yes"
        options["ipopt.warm_start_bound_push"]=1e-6
        options["ipopt.warm_start_slack_bound_push"]=1e-6
        options["ipopt.warm_start_mult_bound_push"]=1e-6
        options["ipopt.mu_init"]=1e-6
        solver = ca.nlpsol("mysolver", Solver, nlp, options)

        solver_in["lbx"]=[-10]*2
        solver_in["ubx"]=[10]*2
        solver_in["lbg"]=[0]
        solver_in["ubg"]=[1]
        solver_in["x0"]=oldsolver_out["x"]
        solver_in["lam_g0"]=oldsolver_out["lam_g"]
        if "bonmin" not in str(Solver): solver_in["lam_x0"] =oldsolver_out["lam_x"]


        solver_out = solver(**solver_in)
        
        if aux_options["codegen"]:
           self.check_codegen(solver,solver_in,digits=codegen_check_digits,**aux_options["codegen"])

  def test_IPOPTrhb2_gen(self):
    self.message("rosenbrock, exact hessian generated, constrained")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    sigma=ca.SX.sym("sigma")
    lambd=ca.SX.sym("lambd")

    for Solver, solver_options, aux_options in solvers:
      if "madnlp" in Solver: continue
      print("test_IPOPTrhb2_gen",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {} #"toldx": 1e-15, "tolgl": 1e-15}).iteritems():
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_out = solver(**solver_in)

      digits = 5

      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(Solver) + str(solver_out["f"][0]) + ":" + str(c_r))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if Solver=="snopt" else 8,str(Solver)+str(6 if Solver=="snopt" else 8))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5 if Solver=="snopt" else 8,str(Solver))
      if Solver not in ["bonmin","sleqp"]: self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))
      self.check_serialize(solver, solver_in)

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,digits=codegen_check_digits,**aux_options["codegen"])


  def test_jacG_empty(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj, 'g':1}

    for Solver, solver_options, aux_options in solvers:
      if "worhp"==Solver:
        continue
      if "sqpmethod"==Solver:
        continue
      if "snopt"==Solver:
        continue
      if "fatrop"==Solver:
        continue
      if "ipmc"==Solver:
        continue  # "Empty sparsity pattern not supported in IPMC C interface"
      if "madnlp"==Solver:
        continue
      print("test_jacG_empty",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[2]
      solver_out = solver(**solver_in)

      digits = 5

      self.checkarray(solver_out["f"],ca.DM([0]),str(Solver),digits=digits)
      self.checkarray(solver_out["x"],ca.DM([1,1]),str(Solver),digits=digits)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0,0]),str(Solver),digits=digits)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],ca.DM([0]),str(Solver),digits=digits)

      if aux_options["codegen"]:
        if "ipopt"==Solver: continue # Empty sparsity pattern not supported in IPOPT C interface
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTrhb2_gen_par(self):
    self.message("rosenbrock, exact hessian generated, constrained, parametric")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    p=ca.SX.sym("p")

    obj = (p-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'p':p, 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTrhb2_gen_par",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_in["p"]=[1]


      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,digits=codegen_check_digits,**aux_options["codegen"])
      solver_out = solver(**solver_in)

      digits = 5

      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5 if Solver=="snopt" else 8,str(Solver))
      if Solver not in ["bonmin","sleqp"]: self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))


  def test_IPOPTrhb_gen(self):
    self.message("rosenbrock, exact hessian generated")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj=(1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj}

    sigma=ca.SX.sym("sigma")

    for Solver, solver_options, aux_options in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      print("test_IPOPTrhb_gen",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      xdig = 6
      
      if Solver=="uno" and "LBFGS" in str(solver_options):
        xdig = 7
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,xdig,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,xdig,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTrhb_gen_xnonfree(self):
    self.message("rosenbrock, exact hessian generated, non-free x")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    obj=(1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'f':obj}

    sigma=ca.SX.sym("sigma")

    for Solver, solver_options, aux_options in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      print("test_IPOPTrhb_gen_xnonfree",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[1,-10]
      solver_in["ubx"]=[1,10]



      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,digits=codegen_check_digits,**aux_options["codegen"])

  def test_IPOPTrhb_gen_par(self):
    self.message("rosenbrock, exact hessian generated, parametric")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    p=ca.SX.sym("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(*[x,y]), 'p':p, 'f':obj}

    sigma=ca.SX.sym("sigma")

    for Solver, solver_options, aux_options in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      print("test_IPOPTrhb_gen_par",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["p"]=1
      solver_out = solver(**solver_in)
      xdig = 6
      if Solver=="uno" and "LBFGS" in str(solver_options):
        xdig = 7
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,xdig,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,xdig,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
        
  @memory_heavy()
  def test_IPOPTnorm(self):
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return c.dot(mx,mx)
    N=10
    x=ca.MX.sym("x",N)
    x0=ca.linspace(0,1,N)
    X0=ca.MX(x0)
    nlp={'x':x, 'f':norm_2(x-X0), 'g':2*x}
    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTnorm",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      # ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():

      solver_in["lbx"]=[-10]*N
      solver_in["ubx"]=[10]*N
      solver_in["lbg"]=[-10]*N
      solver_in["ubg"]=[10]*N
      solver_out = solver(**solver_in)
      print("residuals")
      print(array(solver_out["x"]).squeeze()-x0)
      print("bazmeg", solver_out["f"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.checkarray(array(solver_out["x"]).squeeze(),x0,str(Solver),digits=8)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0]*10),str(Solver),digits=8)
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][1],0,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
  def test_IPOPTnoc(self):
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2}
    for Solver, solver_options, aux_options in solvers:
      if "snopt"==Solver: continue
      print("test_IPOPTnoc",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      # ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTmx(self):
    self.message("trivial IPOPT, using MX")
    x=ca.MX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':2*x}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTmx",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      # ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTc(self):
    self.message("trivial, overconstrained")
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':ca.vertcat(*[x,x,x])}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTc",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10, -10, -10]
      solver_in["ubg"]=[10, 10, 10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,9,str(Solver) )
      self.assertAlmostEqual(solver_out["x"][0],1,5,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTc2(self):
    self.message("trivial2, overconstrained")
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':ca.vertcat(*[x,x,x+x])}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTc2",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10, -10, -10]
      solver_in["ubg"]=[10, 10, 10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
  def test_IPOPTcmx(self):
    self.message("trivial , overconstrained, using MX")
    x=ca.MX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':ca.vertcat(*[2*x,3*x,4*x])}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTcmx",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10,-10,-10]
      solver_in["ubg"]=[10,10,10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,9,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    nlp={'x':ca.vertcat(*[x,y]), 'f':0, 'g':ca.vertcat(*[x-y,x])}
    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTdeg",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10, -10]
      solver_in["ubx"]=[10, 10]
      solver_in["lbg"]=[0, 3]
      solver_in["ubg"]=[0, 3]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],solver_out["x"][1],4 if "sqic" in str(solver_options) else 10,"IPOPT")

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_IPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    nlp={'x':ca.vertcat(*[x,y]), 'f':0, 'g':ca.vertcat(*[x-y,x,x+y])}

    for Solver, solver_options, aux_options in solvers:
      print("test_IPOPTdegc",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10, -10]
      solver_in["ubx"]=[10, 10]
      solver_in["lbg"]=[0, 3 , -10]
      solver_in["ubg"]=[0, 3, 10]
      solver_out = solver(**solver_in)
      # todo: catch error when set([0, 3 , 5]) two times
      self.assertAlmostEqual(solver_out["x"][0],solver_out["x"][1],4 if "sqic" in str(solver_options) else 10,"IPOPT")

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_XfreeChange(self):
    self.message("Change in X settings")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_XfreeChange",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      if "qrqp" in str(solver_options): solver_in["x0"]=[0.6,1]
      solver_in["lbx"]=[-10,-10]
      solver_in["ubx"]=[10,10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      solver_in["lbx"]=[-10,1]
      solver_in["ubx"]=[10,1]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]



      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["f"][0],0,9,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))
      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_activeLBX(self):
    self.message("active LBX")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_activeLBX",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      if "qrqp" in str(solver_options): solver_in["x0"]=[0.5,1]
      solver_in["lbx"]=[-10,1.2]
      solver_in["ubx"]=[10,2]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      if float(solver_out["x"][0])<0: # JOEL: There appears to be two local minima
        self.assertAlmostEqual(solver_out["f"][0],4.3817250416084308,6,str(Solver))
        self.assertAlmostEqual(solver_out["x"][0],-1.0910624688699295,6,str(Solver))
        self.assertAlmostEqual(solver_out["x"][1],1.2,5,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if "stabilizedsqp"==Solver else 8,str(Solver)+str(solver_options))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],-1.9165378046901287,4,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,8,str(Solver))
      else:
        self.assertAlmostEqual(solver_out["f"][0],9.0908263002590e-3,6,str(Solver))
        self.assertAlmostEqual(solver_out["x"][0],1.0952466252248,6,str(Solver))
        self.assertAlmostEqual(solver_out["x"][1],1.2,5,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if "stabilizedsqp"==Solver else 8,str(Solver)+str(solver_options))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],-8.6963632695079e-2,4,str(Solver))
        if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
  def test_activeLBG(self):
    self.message("active LBG")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_activeLBG",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      if "qrqp" in str(solver_options): solver_in["x0"]=[1.5,1]
      solver_in["lbx"]=[-10,-10]
      solver_in["ubx"]=[10,10]
      solver_in["lbg"]=[2.2]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],4.252906468284e-3,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1.065181061847138,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1.1348189166291160,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,4,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],-4.1644422845712e-2,3,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_activeUBG(self):
    self.message("active UBG")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_activeUBG",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      if "qrqp" in str(solver_options): solver_in["x0"]=[1.5,1]
      solver_in["lbx"]=[-10,-10]
      solver_in["ubx"]=[10,10]
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1.8]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],4.64801220074552e-3,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],9.318651964592811e-1,5,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],8.68134821123689e-1,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,4,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],4.75846495145007e-2,5,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_activeUBX(self):
    self.message("active UBX")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, aux_options in solvers:
      print("test_activeUBX",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0,1]
      solver_in["lbx"]=[-10,0]
      solver_in["ubx"]=[10,0.9]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],2.626109721583e-3,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],9.4882542279172277e-01,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],0.9,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],5.39346608659e-2,4,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,8,str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
  @memory_heavy()
  def test_QP(self):
    self.message("QP")

    N = 50

    x = ca.SX.sym("x",N)
    x0 = ca.DM(list(range(N)))
    H = ca.diag(list(range(1,N+1)))
    obj = 0.5*ca.mtimes([(x-x0).T,H,(x-x0)])

    nlp = {'x':x, 'f':obj}
    for Solver, solver_options, aux_options in solvers:
      if "snopt"==Solver: continue
      if "fatrop"==Solver: continue
      if Solver=="sqpmethod" and "limited-memory" in str(solver_options): continue
      print("test_QP",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=-1000
      solver_in["ubx"]=1000
      solver_out = solver(**solver_in)
      self.checkarray(solver_out["x"],x0,str(Solver),digits=2)
      self.assertAlmostEqual(solver_out["f"][0],0,3,str(Solver))
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM.zeros(N,1),str(Solver),digits=4)

  def test_QP2(self):
    H = ca.DM([[1,-1],[-1,2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([[1, 1],[-1, 2],[2, 1]])

    LBA = ca.DM([-inf]*3)
    UBA = ca.DM([2, 2, 3])

    LBX = ca.DM([0.5,0])
    UBX = ca.DM([0.5,inf])

    x=ca.SX.sym("x",2)
    nlp={'x':x, 'f':0.5*ca.mtimes([x.T,H,x])+G.T @ x, 'g':A @ x}

    for Solver, solver_options, aux_options in solvers:
      print("test_QP2",Solver,solver_options)
      options = dict(solver_options)
      if "madnlp" in str(Solver): continue
      if "ipopt" in str(Solver):
        options["ipopt.fixed_variable_treatment"] = "make_constraint"
      solver = ca.nlpsol("mysolver", Solver, nlp, options)
      #{"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
      solver_in = {}

      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA
      if 'sqic' in str(solver_options):
        continue
      if 'daqp' in str(solver_options):
        continue

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1.25,6,str(Solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],4.75,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],ca.DM([0,2,0]),str(Solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-7.4375,6,str(Solver))
      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

      solver = ca.nlpsol("mysolver", Solver, nlp, options)
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1.25,6,str(Solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],4.75,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],ca.DM([0,2,0]),str(Solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-7.4375,6,str(Solver))
      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_QP2_unconvex(self):
    H = ca.DM([[1,-1],[-1,-2]])
    G = ca.DM([-2,-6])
    A =  ca.DM([[1, 1],[-1, 2],[2, 1]])

    LBA = ca.DM([-inf]*3)
    UBA = ca.DM([2, 2, 3])

    LBX = ca.DM([0]*2)
    UBX = ca.DM([inf]*2)

    x=ca.SX.sym("x",2)
    nlp={'x':x, 'f':0.5*ca.mtimes([x.T,H,x])+G.T @ x, 'g':A @ x}

    for Solver, solver_options, aux_options in solvers:
      print("test_QP2_unconvex",Solver,solver_options)
      options = dict(solver_options)
      if "ipopt" in str(Solver):
        options["ipopt.fixed_variable_treatment"] = "make_constraint"
      solver = ca.nlpsol("mysolver", Solver, nlp, options)
      solver_in = {}
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      if "qrqp" in str(solver_options):
        solver_in["x0"]= ca.DM([1,1])
        solver_in["lam_g0"]= ca.DM([1,1,0])
      if 'daqp' in str(solver_options):
        continue
      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,6,str(solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(solver))

      if "bonmin" not in str(Solver) and "worhp" not in str(Solver): self.checkarray(solver_out["lam_g"],ca.DM([4+8.0/9,20.0/9,0]),str(solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-10-16.0/9,6,str(solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])
      solver = ca.nlpsol("mysolver", Solver, nlp, options)

      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,6,str(solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(solver))

      if "bonmin" not in str(Solver) and "worhp" not in str(Solver): self.checkarray(solver_out["lam_g"],ca.DM([4+8.0/9,20.0/9,0]),str(solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-10-16.0/9,6,str(solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_bug(self):
    x = ca.MX.sym("x", 3)
    y = ca.MX.sym("y", 2)
    f = ca.Function("f", [x, y], [1.])

    aa = ca.MX.sym("aa", 5)
    a = aa[:3]
    b = aa[3:]
    f_call = f(a, b)
    nlp = {'x':aa, 'f':f_call}
    for Solver, solver_options, aux_options in solvers:
      if "worhp" in Solver: continue
      if "snopt"==Solver: continue
      print("test_bug",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

  def test_missing_symbols(self):
    x = ca.MX.sym("x")
    p = ca.MX.sym("p")

    for Solver, solver_options, aux_options in solvers:
      print("test_missing_symbols",Solver,solver_options)
      with self.assertInException("[p] are free"):
        solver = ca.nlpsol("solver",Solver,{"x":x,"f":(x-p)**2}, solver_options)

  @requires_nlpsol("ipopt")
  def test_no_success(self):

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    f = (1-x)**2+100*(y-x**2)**2
    for Solver, solver_options, aux_options in solvers:
      print("test_no_success",Solver,solver_options)
      solver = ca.nlpsol("solver","ipopt",{'x':ca.vertcat(x,y), 'f':f,'g':ca.vertcat(x+1,x-2)})
      solver(x0=0,lbg=0,ubg=0)
      self.assertFalse(solver.stats()["success"])

      solver = ca.nlpsol("solver","ipopt",{'x':ca.vertcat(x,y), 'f':f,'g':ca.vertcat(x+1,x-2)},{"error_on_fail":True})
      with self.assertInException("process"):
        solver(x0=0,lbg=0,ubg=0)

  
  @requires_nlpsol("ipopt")
  def test_postpone_expand(self):
    x = ca.MX.sym("x")
    p = ca.MX.sym("p")
    
    J = ca.Function("jac_f", [x, ca.MX(1,1)], [-x], ['x', 'out_r'], ['jac_r_x'])
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":f(x-3)},{"ipopt.max_iter": 5})
    solver()
    self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_LIMITED")
    self.assertTrue(solver.get_function('nlp_jac_g').is_a("MXFunction"))

    J = ca.Function("jac_f", [x, ca.MX(1,1)], [-x], ['x', 'out_r'], ['jac_r_x'])
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":f(x-3)},{"ipopt.max_iter": 5, "expand":True})
    solver()
    self.assertTrue(solver.stats()["return_status"]=="Solve_Succeeded")
    self.assertTrue(solver.get_function('nlp_jac_g').is_a("SXFunction"))
    


    J = ca.Function("jac_f", [x, ca.MX(1,1)], [-x], ['x', 'out_r'], ['jac_r_x'])
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":f(x-3)},{"ipopt.max_iter": 5, "postpone_expand":True})
    solver()
    self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_LIMITED")
    self.assertTrue(solver.get_function('nlp_jac_g').is_a("MXFunction"))

    J = ca.Function("jac_f", [x, ca.MX(1,1)], [-x], ['x', 'out_r'], ['jac_r_x'])
    f = ca.Function('f', [x], [x**2], ['x'], ['r'], dict(custom_jacobian = J, jac_penalty = 0))
    
    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":f(x-3)},{"ipopt.max_iter": 5, "expand":True, "postpone_expand":True})
    solver()
    self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_LIMITED")
    self.assertTrue(solver.get_function('nlp_jac_g').is_a("SXFunction"))
    

  @requires_nlpsol("ipopt")
  def test_iteration_Callback(self):

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    f = (1-x)**2+100*(y-x**2)**2
    nlp={'x':ca.vertcat(x,y), 'f':f,'g':x+y}
    fcn = ca.Function('f', [x, y], [f])

    class MyCallback(ca.Callback):
      def __init__(self,nx, ng, np):
        ca.Callback.__init__(self)
        self.foo = []

        self.nx = nx
        self.ng = ng
        self.np = np
        self.construct("mycallback", {})

      def get_n_in(self): return ca.nlpsol_n_out()
      def get_n_out(self): return 1


      def get_sparsity_in(self, i):
        n = ca.nlpsol_out(i)
        if n=='f':
          return ca.Sparsity. scalar()
        elif n in ('x', 'lam_x'):
          return ca.Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
          return ca.Sparsity.dense(self.ng)
        else:
          return ca.Sparsity(0,0)
      def eval(self, arg):
        self.foo.append(arg)
        return [0]

    mycallback = MyCallback(2,1,0)
    opts = {}
    opts['iteration_callback'] = mycallback
    opts['ipopt.tol'] = 1e-8
    opts['ipopt.max_iter'] = 50
    solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
    sol = solver(lbx=-10, ubx=10, lbg=-10, ubg=10)
    self.assertEqual(len(mycallback.foo),solver.stats()["iter_count"]+1)

    class MyCallback(ca.Callback):
      def __init__(self,nx, ng, np):
        ca.Callback.__init__(self)
        self.foo = []

        self.nx = nx
        self.ng = ng
        self.np = np
        self.construct("mycallback", {})

      def get_n_in(self): return ca.nlpsol_n_out()
      def get_n_out(self): return 1


      def get_sparsity_in(self, i):
        n = ca.nlpsol_out(i)
        if n=='f':
          return ca.Sparsity. scalar()
        elif n in ('x', 'lam_x'):
          return ca.Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
          return ca.Sparsity.dense(self.nx)
        else:
          return ca.Sparsity(0,0)
      def eval(self, arg):
        self.foo.append(arg)
        return [0]

    mycallback = MyCallback(2,1,0)
    opts = {}
    opts['iteration_callback'] = mycallback
    opts['ipopt.tol'] = 1e-8
    opts['ipopt.max_iter'] = 50

    try:
      solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
    except Exception as e:
      self.assertTrue("Callback function input size mismatch" in str(e))
  def test_pathological(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+y**2}

    for Solver, solver_options, aux_options in solvers:
      if "snopt"==Solver: continue
      if "worhp"==Solver or "stabilizedsqp"==Solver : continue
      print("test_pathological",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,-1]
      solver_in["ubx"]=[10,2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],ca.DM([0]),digits=7)
      self.checkarray(solver_out["x"],ca.DM([1,0]),digits=7,failmessage=str(Solver))
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0,-0]),digits=7,failmessage=str(Solver))

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_pathological2(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2+y}

    for Solver, solver_options, aux_options in solvers:
      if "snopt"==Solver: continue
      print("test_pathological2",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,0]
      solver_in["ubx"]=[10,2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],ca.DM([0]),digits=7)
      self.checkarray(solver_out["x"],ca.DM([1,0]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0,-1]),digits=7)

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_pathological3(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    nlp={'x':ca.vertcat(*[x,y]), 'f':(1-x)**2, 'g':x+y}

    for Solver, solver_options, aux_options in solvers:
      if "worhp"==Solver: continue
      print("test_pathological3",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,0]
      solver_in["ubx"]=[10,2]
      solver_in["lbg"]=[2]
      solver_in["ubg"]=[2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],ca.DM([0]),digits=7)
      self.checkarray(solver_out["x"],ca.DM([1,1]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0,0]),digits=7)

      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_pathological4(self):
    x=ca.SX.sym("x")
    nlp={'x':x, 'f':x*x}

    for Solver, solver_options, aux_options in solvers:
      if "snopt"==Solver: continue
      if "worhp"==Solver: continue
      if "madnlp"==Solver: continue
      print("test_pathological4",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[0]
      solver_in["lbx"]=[0]
      solver_in["ubx"]=[0]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],ca.DM([0]),digits=7)
      self.checkarray(solver_out["x"],ca.DM([0]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],ca.DM([0]),digits=7)
      if aux_options["codegen"]:
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  def test_nlp_sensitivity(self):

    x = ca.MX.sym("x")
    p = ca.MX.sym("p")

    nlp = {"x":x,"p":p,"f":(ca.sin(x)-p**2)**2,"g":x}

    for Solver, solver_options, aux_options in solvers:
      if "madnlp" in str(Solver): continue
      if "fatrop" in str(Solver): continue
      if "ipmc" in str(Solver): continue
      if "ipopt" in str(solver_options): continue
      if "snopt" in str(solver_options): continue
      if Solver in ["alpaqa"]: continue
      if "knitro" in str(Solver):
        solver_options = copy.deepcopy(solver_options)
        solver_options["knitro"]["algorithm"] = 4 # sqp
      print("test_nlp_sensitivity",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      print("foo",Solver,solver(lbg=0,x0=0,p=0.5))

      z = solver(p=p,x0=x,lbg=0)["x"]

      z2 = ca.asin(p**2)

      f = ca.Function('f',[x,p],[z,ca.jacobian(z,p)])
      f2 = ca.Function('f',[x,p],[z2,ca.jacobian(z2,p)])

      self.checkfunction_light(f,f2,[0,0.5],digits=6)

  @requires_conic("qrqp")
  def test_regularize_sqpmethod(self):

    # Test problem that is indefinite in direction of the constraint Jacobian
    x = ca.MX.sym("x",2)
    f = 0.5*ca.bilin(ca.DM([[1,0],[0,-2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    solver = ca.nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True}})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats = solver.stats()
    self.assertTrue(stats["iter_count"]==1)
    self.assertTrue("H:\n[[1, 0], \n [0, -2]]" in result[0])
    self.checkarray(res["x"],ca.DM([0,2]),digits=6)

    solver = ca.nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize","convexify_margin":0})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],ca.DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==2)
    self.assertTrue("H:\n[[3, 0], \n [0, 0]]" in result[0])

    solver = ca.nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize","convexify_margin":1e-4})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],ca.DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==2)
    self.assertTrue("H:\n[[3.0001, 0], \n [0, 0.0001]]" in result[0])

    x = ca.MX.sym("x",2)
    f = 0.5*ca.bilin(ca.DM([[1,0],[0,2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    solver = ca.nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True}})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats = solver.stats()
    self.assertTrue(stats["iter_count"]==1)
    print(result[0])
    self.assertTrue("H:\n[[1, 0], \n [0, 2]]" in result[0])
    self.checkarray(res["x"],ca.DM([0,2]),digits=6)

    solver = ca.nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize"})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],ca.DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==1)
    self.assertTrue("H:\n[[1, 0], \n [0, 2]]" in result[0])

  def test_infeasible(self):
    x = ca.MX.sym("x")

    nlp = {"f":x**2,"x":x,"g":ca.sin(x)}

    for Solver, solver_options, aux_options in solvers:
        myoptions = dict(solver_options)
        myoptions["error_on_fail"] = True
        if Solver=="sqpmethod": continue
        solver = ca.nlpsol("solver",Solver,nlp,myoptions)
        
        if aux_options["codegen"] and args.run_slow:
          cg = self.check_codegen(solver,dict(lbg=-10,ubg=10),**aux_options["codegen"])
          F = cg["F"]
          with self.assertRaises(RuntimeError):
            F(lbg=5,ubg=10)
          
          if os.name != 'nt':
              # uno's interior-point ("ipopt") preset reports an infeasible problem
              # as an algorithmic error (line-search failure) -> SOLVER_RET_EXCEPTION,
              # whereas the SQP presets report infeasible-stationary. Both are valid
              # "solver failed on infeasible problem" outcomes (vm and codegen agree).
              self.check_codegen(solver,dict(lbg=5,ubg=10),main=True,main_return_code=[ca.SOLVER_RET_UNKNOWN, ca.SOLVER_RET_INFEASIBLE, ca.SOLVER_RET_EXCEPTION],**aux_options["codegen"])

  def test_indefinite(self):

    # Test problem that is indefinite in direction of the constraint Jacobian
    x = ca.MX.sym("x",2)
    f = 0.5*ca.bilin(ca.DM([[1,0],[0,-2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    for Solver, solver_options, aux_options in solvers:
      solver_in = {"lbg": 2, "ubg": 2}
      print("test_indefinite",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      out = solver(**solver_in)
      self.checkarray(out["x"],ca.DM([0,2]),digits=6)
      if "bonmin" not in str(Solver): self.checkarray(out["lam_g"],ca.DM([4]),digits=6)

      if aux_options["codegen"]:
        solver.generate('f.c',{"main":True})
        solver.generate_in("in.dat",solver.convert_in(solver_in))
        print(solver_in)
        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  @requires_nlpsol("sqpmethod")
  def test_gauss_newton_sqpmethod(self):
    x = ca.SX.sym("x",3)

    F = ca.sin(x) - ca.vertcat(1,2,3)*ca.vertcat(x[2],0,0)
    J = ca.jacobian(F,x)
    f = 0.5*ca.dot(F,F)
    p = ca.SX.sym("x",0,1)
    lam_f = ca.SX.sym("x")
    lam_g = ca.SX.sym("x",0,1)
    GN = ca.Function('GN',[x,p,lam_f,lam_g],[lam_f*ca.triu(J.T @ J)])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    with self.assertInException("Hessian must be symmetric"):
      solver = ca.nlpsol("solver","sqpmethod",nlp,options)

    # A 2-norrm problem ...
    F = ca.sin(x-ca.vertcat(0.22,0.72,0.2)) - ca.vertcat(0.1,0.5,0.99)
    J = ca.jacobian(F,x)
    f = 0.5*ca.dot(F,F)
    H = ca.Function("H",[x],[ca.hessian(f,x)[0]])
    H = H(0)
    # with an indefinite Hessian at x0
    self.assertTrue(np.any(np.linalg.eigh(H)[0]<0))

    # Solve with Gauss-Newton -> 6 iterations
    GN = ca.Function('GN',[x,p,lam_f,lam_g],[lam_f*(J.T @ J)])
    options = {"convexify_strategy":"regularize","qpsol":"qrqp","hess_lag": GN}
    nlp = {"x":x,"f":f}
    solver = ca.nlpsol("solver","sqpmethod",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==6)

    # Solve with exact Hessian + regularization -> 9 iterations
    options = {"convexify_strategy":"regularize","qpsol":"qrqp"}
    nlp = {"x":x,"f":f}
    solver = ca.nlpsol("solver","sqpmethod",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==9)

  @requires_nlpsol("ipopt")
  def test_gauss_newton_ipopt(self):
    x = ca.SX.sym("x",3)

    F = ca.sin(x) - ca.vertcat(1,2,3)*ca.vertcat(x[2],0,0)
    J = ca.jacobian(F,x)
    f = 0.5*ca.dot(F,F)
    p = ca.SX.sym("x",0,1)
    lam_f = ca.SX.sym("x")
    lam_g = ca.SX.sym("x",0,1)
    GN = ca.Function('GN',[x,p,lam_f,lam_g],[lam_f*(J.T @ J)])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    with self.assertInException("Hessian must be upper triangular"):
      solver = ca.nlpsol("solver","ipopt",nlp,options)

    # A 2-norrm problem ...
    F = ca.sin(x-ca.vertcat(0.22,0.72,0.2)) - ca.vertcat(0.1,0.5,0.99)
    J = ca.jacobian(F,x)
    f = 0.5*ca.dot(F,F)
    H = ca.Function("H",[x],[ca.hessian(f,x)[0]])
    H = H(0)
    # with an indefinite Hessian at x0
    self.assertTrue(np.any(np.linalg.eigh(H)[0]<0))

    # Solve with Gauss-Newton -> 6 iterations
    GN = ca.Function('GN',[x,p,lam_f,lam_g],[lam_f*ca.triu(J.T @ J)])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    solver = ca.nlpsol("solver","ipopt",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==6)

    # Solve with exact Hessian + regularization -> 9 iterations
    nlp = {"x":x,"f":f}
    solver = ca.nlpsol("solver","ipopt",nlp)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==9)

  def test_cvx_sqpmethod(self):

    eps = 1e-2

    H = ca.diagcat(ca.DM([[2,3],[3,0]]),ca.DM([[7,1,0],[1,0,6],[0,6,3]]))
    p = [0,3,2,1,4]
    Hp = H[p,p]

    eig = np.linalg.eig

    def reflect(A,eps):
      [d,V] = np.linalg.eigh(A)  
      d = abs(d)
      D = ca.diag(d)
      return np.dot(V,np.dot(D,V.T))

    def clip(A,eps):
      [d,V] = np.linalg.eigh(A)
      d[d<eps] = eps
      D = ca.diag(d)  # pyright: ignore[reportCallIssue,reportArgumentType]
      return np.dot(V,np.dot(D,V.T))


    for H, Hcvx, opts, includes_init, excludes_init in [
      (ca.sparsify(ca.DM([[0, 0],[0, 2]])), ca.DM([[eps, 0],[0, 2+eps]]), {"convexify_strategy": "regularize","convexify_margin":eps},[],[]),
      (ca.diagcat(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]]),ca.DM([[5, 3],[3, 6]])), ca.diagcat(ca.DM([[2, 0,0],[0, 1,0],[0,0,8]]),ca.DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 2 blocks with maximum size 3"],[]),
      (ca.diagcat(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]]),ca.DM([[5, 3],[3, 6]])), ca.diagcat(ca.DM([[2, 0,0],[0, eps,0],[0,0,8]]),ca.DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      (ca.DM([[2, 0,0],[0, eps/2,0],[0,0,8]]), ca.DM([[2, 0,0],[0, eps,0],[0,0,8]]), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (ca.diagcat(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]]),ca.DM([[5, 8],[8, 6]])), ca.diagcat(ca.DM([[2, 0,0],[0, 1,0],[0,0,8]]),reflect(ca.DM([[5, 8],[8, 6]]),eps)), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (ca.diagcat(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]]),ca.DM([[5, 8],[8, 6]])), ca.diagcat(ca.DM([[2, 0,0],[0, eps,0],[0,0,8]]),clip(ca.DM([[5, 8],[8, 6]]),eps)), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      (ca.DM([[2, 0,0],[0, eps/2,0],[0,0,8]]), ca.DM([[2, 0,0],[0, eps,0],[0,0,8]]), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (ca.diagcat(ca.sparsify(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]])),ca.DM([[5, 3],[3, 6]])), ca.diagcat(ca.DM([[2, 0,0],[0, 1,0],[0,0,8]]),ca.DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 4 blocks with maximum size 2"],[]),
      (ca.diagcat(ca.sparsify(ca.DM([[2, 0,0],[0, -1,0],[0,0,8]])),ca.DM([[5, 3],[3, 5]])), ca.diagcat(ca.DM([[2, 0,0],[0, 1,0],[0,0,8]]),ca.DM([[5, 3],[3, 5]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 4 blocks with maximum size 2"],[]),
      (Hp, reflect(Hp,eps), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (Hp, clip(Hp,eps), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      ]:

      n = H.shape[0]
      x = ca.MX.sym("H",n)
      nlp = {"x": x, "f": 0.5*ca.bilin(H,x,x)}
      options = {"max_iter":1,"qpsol":"qrqp","verbose":True}
      options.update(opts)
      with self.assertOutput(includes_init,excludes_init): 
        solver = ca.nlpsol("solver","sqpmethod",nlp,options)

      x0 = ca.DM.ones(n,1)

      res = solver(x0=x0)

      self.checkarray(x0-np.linalg.solve(np.array(Hcvx),np.array(H @ x0)),res["x"])

      self.check_serialize(solver,{"x0":x0})

      self.check_codegen(solver,{"x0":x0},std="c99",digits=codegen_check_digits)


  def test_simple_bounds_detect(self):

    x = ca.SX.sym("x",5)
    p = ca.SX.sym("p",1)


    g = [
      (1.1,  x[0]*x[1], 2),
      (-inf,  x[4], 2), # 4 H
      (-10,  x[0], 10),
      (-5,  x[0], 2), # 0 H
      (-4,  x[0], 4), # 0 L
      (1.1,  x[4]*x[1], 2),
      (0,  x[4], inf), # 4 L
      (7,  x[2], 7), # 2 LH
      (-4,  x[2], 40),
      (9,  x[1], 9), # 1 LH
      (-4,  x[0], 4), # 0 L
      (-4,  x[1], 9)] # 1 H

    [lbg,g,ubg]= zip(*g)

    lbg = ca.vcat(lbg)
    ubg = ca.vcat(ubg)
    g = ca.vcat(g)

    [gi,lbx,ubx,lam_f,lam_b]=ca.detect_simple_bounds(x,p,g,lbg,ubg)

    self.checkarray(ca.DM(lbx).T,[-4,9,7,-inf,0])
    self.checkarray(ca.DM(ubx).T,[2,9,7,inf,2])

    def round_trip_f(arg):
      return lam_b(*(lam_f(arg,0)+(0,)))

    def round_trip_b(arg):
      return lam_f(lam_b(*(arg+(0,))),0)


    G = np.array([[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]])

    for i in range(12):
      a = ca.DM.zeros(12,1)
      a[i] = 2
      b = round_trip_f(a)
      c = round_trip_f(round_trip_f(a))
      self.checkarray(b,G[i,:])
      self.checkarray(b,c)

    G = np.array([[-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0],
                  [0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0],
                  [0, 0, 0, 0,-1, 0, 0, 0, 0, 0,-1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


    for i in range(12):
      a = ca.DM.zeros(12,1)
      a[i] = -2
      b = round_trip_f(a)
      c = round_trip_f(round_trip_f(a))
      self.checkarray(b,G[i,:])
      self.checkarray(b,c)

    G = np.array([[-2, 0, 0, 0, 0],
                  [0, -2, 0, 0, 0],
                  [0, 0, -2, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, -2]])


    for i in range(5):
      a = ca.DM.zeros(5,1)
      a[i] = -2
      b = round_trip_b(([3,7],a))
      c = round_trip_b(round_trip_b(([3,7],a)))
      self.checkarray(b[0].T,[3,7])
      self.checkarray(b[1],G[i,:])
      self.checkarray(b[1],c[1])

    G = np.array([[2, 0, 0, 0, 0],
                  [0, 2, 0, 0, 0],
                  [0, 0, 2, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 2]])


    for i in range(5):
      a = ca.DM.zeros(5,1)
      a[i] = 2
      b = round_trip_b(([3,7],a))
      c = round_trip_b(round_trip_b(([3,7],a)))
      self.checkarray(b[0].T,[3,7])
      self.checkarray(b[1],G[i,:])
      self.checkarray(b[1],c[1])
      
  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_simple_bounds_detect_solvers(self):
  
    for Solver, solver_options, aux_options in solvers:
        print(Solver,solver_options)
        #if Solver=="bonmin": continue
        if Solver=="sleqp": continue
        if Solver=="madnlp": continue

        # detect_simple_bounds turns "3 <= x[2] <= 3" and "2 <= x[1] <= 2"
        # below into FIXED variables, and ipopt's default
        # fixed_variable_treatment ("make_parameter") drops fixed variables
        # from the NLP and reports a zero multiplier for them.  The reference
        # solve, where those rows stay general constraints, reports the real
        # ones (-6, 7.6, 34, -46 depending on the case), so the two disagree
        # in lam_g by up to 46 for reasons that have nothing to do with
        # detect_simple_bounds.  Ask ipopt to keep them.
        if Solver == "ipopt":
            solver_options = dict(solver_options)
            solver_options["ipopt"] = dict(solver_options.get("ipopt", {}),
                                           fixed_variable_treatment="make_constraint")

      
        x = ca.MX.sym("x",5)
        
        # different a, b


        g = [
          (1.1,  x[0]*x[1], 2),
          (-11,  x[4], 2), # 4 H
          (-10,  x[0], 10),
          (-5,  x[0], 2), # 0 H
          (-4,  x[0], 4), # 0 L
          (1.1,  x[4]*x[1], 2),
          (0,  x[4], inf), # 4 L
          (3,  x[2], 3), # 2 LH
          (-4,  x[2], 40),
          (2,  x[1], 2), # 1 LH
          (-4,  x[0], 4), # 0 L
          (-4,  x[1], 9)] # 1 H

        [lbg,g,ubg]= zip(*g)

        lbg = ca.vcat(lbg)
        ubg = ca.vcat(ubg)
        g = ca.vcat(g)
        
        I = ca.densify(ca.DM.eye(5))

        for lbx,ubx in [(-ca.DM.inf(5),ca.DM.inf(5))]:#,(-3*DM.ones(5),3*DM.ones(5)),(-4*DM.ones(5),4*DM.ones(5))]:
        
            for i in range(3):#range(5):
                for sign in [-1,1]:
                    x_target = I[:,i]*sign*20
                    
                    nlp = {"x":x,"g":g,"f":ca.sumsqr(x-x_target)}
                    
                    solver_ref_ipopt = ca.nlpsol("mysolver", "ipopt", nlp)
                    
                    solver_ref = ca.nlpsol("mysolver", Solver, nlp, solver_options)
                    
                    
                    
                    my_solver_options = dict(solver_options)
                    my_solver_options["detect_simple_bounds"] = True
                    
                    solver = ca.nlpsol("mysolver", Solver, nlp, my_solver_options)
                    solver_in = dict(lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
                    
                    sol_ref_ipopt = solver_ref_ipopt(**solver_in)
                    
                    solver_in["x0"] = sol_ref_ipopt["x"]*1.2
                    
                    digits = 6
                    if "daqp" in str(solver_options):
                        #digits = 4
                        continue
                    if "worhp" in str(solver_options):
                        digits = 4
                        
                    solver_ref(**solver_in)
                    
                    print("stats",solver_ref.stats())

                    if Solver=="bonmin":
                        solver_ref_out = solver_ref(**solver_in)
                        solver_out = solver(**solver_in)
                        for output in ["x","f"]:
                            self.checkarray(solver_out[output],solver_ref_out[output],digits=digits,failmessage=str(Solver)+str(solver_options))
                    else:
                        self.checkfunction_light(solver, solver_ref, inputs=solver_in,digits=digits,failmessage=str(Solver)+str(solver_options))
                    
                    if aux_options["codegen"]:
                        self.check_codegen(solver,solver_in,**aux_options["codegen"])

  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_simple_bounds_sign_issue(self):
  
    for X in [ca.SX,ca.MX]:
        x = X.sym("x")
        
        for f in [x,-x]:
        
            
            for lbg,g,ubg in [
                (-2,x,3),
                (2,x,2),
                (-2,0.3*x,3),
                (2,0.7*x,2),
                (-2,-x,3),
                (-2,-0.7*x,3),
                (-2,0.3*x,3),
                ]:
                
              # Both options are needed to make the two formulations
              # comparable at the default digits=9, and neither has anything
              # to do with detect_simple_bounds itself:
              #  * bound_relax_factor: ipopt widens a simple BOUND by
              #    1e-8*max(1,|bound|) before starting the barrier and never
              #    undoes it, but leaves a general constraint row alone, so
              #    the detected solver lands ~2e-8 outside the bound the
              #    reference sits exactly on;
              #  * fixed_variable_treatment: for the (2, x, 2) rows detection
              #    produces lbx == ubx, and ipopt's default drops such
              #    variables and reports lam = 0 for them (lam_g off by 1.0
              #    and 1.43 on those two cases).
              ipopt_opts = {"ipopt": {"bound_relax_factor": 0,
                                      "fixed_variable_treatment": "make_constraint"}}
              solver_ref = ca.nlpsol("solver","ipopt",{"x":x,"f":f,"g":g},dict(ipopt_opts))
              solver = ca.nlpsol("solver","ipopt",{"x":x,"f":f,"g":g},dict(ipopt_opts,detect_simple_bounds=True))
              
              self.checkfunction_light(solver,solver_ref,inputs=solver.convert_in(dict(lbg=lbg,ubg=ubg)))

  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_simple_bounds_equality(self):
  
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.vertcat(x,y)
    f = x
    
    g = ca.vertcat(x**2+y**2,x,x*y)
    lbg = ca.vertcat(1,-2,-5)
    ubg = ca.vertcat(1,2,5)
    
    x0 = ca.vertcat(0.2,0.6)
    
    
    solver_ref = ca.nlpsol("solver","ipopt",{"x":z,"f":f,"g":g},{"equality": [True,False,False]})
    solver = ca.nlpsol("solver","ipopt",{"x":z,"f":f,"g":g},{"detect_simple_bounds":True, "equality":[True,False,False]})
    
    self.checkfunction_light(solver,solver_ref,inputs=solver.convert_in(dict(x0=x0,lbg=lbg,ubg=ubg)))
  
  @requires_nlpsol("ipopt")
  def test_issue_3407(self):
    opti = ca.Opti()
    x = opti.variable()

    opti.minimize(x**2)

    f = ca.Function('f',[x],[3*x,x.attachAssert(x>1)**2],{"never_inline":True})

    y,z = f(x)

    opti.subject_to(y<=10)
    opti.subject_to(z==3)

    opti.set_initial(x, 3)

    opti.solver("ipopt",{"detect_simple_bounds":True})

    sol = opti.solve()

  @requires_nlpsol("ipopt")
  def test_issue_4235(self):
    opti = ca.Opti()
    x = opti.variable()

    opti.minimize(x**2)

    f = ca.Function('f',[x],[3*x,x.attachAssert(x>1000)**2])

    y,z = f(2*x)

    opti.subject_to(y<=10)
    opti.subject_to(z==3)

    opti.set_initial(x, 3)

    opti.solver("ipopt",{"detect_simple_bounds":True})

    with self.assertRaises(Exception):
        sol = opti.solve()

    opti = ca.Opti()
    x = opti.variable()

    opti.minimize(x**2)

    f = ca.Function('f',[x],[3*x,x.attachAssert(x>1000)**2])

    y,z = f(2*x)

    opti.subject_to(y<=10)
    opti.subject_to(z==3)

    opti.set_initial(x, 3)

    # Expand will have the side effect of removing the assertion
    opti.solver("ipopt",{"expand":True,"detect_simple_bounds":True})

    sol = opti.solve()
        
  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_simple_bounds_detect2(self):
    x = ca.MX.sym("x",5)
    p = ca.MX.sym("p",5)
    
    # different a, b


    g = [
      (1.1,  x[0]*x[1], 2),
      (-11,  x[4], 2), # 4 H
      (-10,  x[0], 10),
      (-5,  x[0], 2), # 0 H
      (-4,  x[0], 4), # 0 L
      (1.1,  x[4]*x[1], 2),
      (0,  x[4], inf), # 4 L
      (3,  x[2], 3), # 2 LH
      (-4,  x[2], 40),
      (2,  x[1], 2), # 1 LH
      (-4,  x[0], 4), # 0 L
      (-4,  x[1], 9)] # 1 H

    [lbg,g,ubg]= zip(*g)

    lbg = ca.vcat(lbg)
    ubg = ca.vcat(ubg)
    g = ca.vcat(g)
    
    def merge(a,b):
        a = dict(a)
        a.update(b)
        return a
        

    nlp = {"x":x,"p":p,"g":g,"f":ca.sumsqr(x-p)}
    
    print_level = 0
    
    with self.assertOutput(["Total number of inequality constraints...............:        2"],[]):
        ca.nlpsol("solver","ipopt",nlp,{"detect_simple_bounds": True})(x0=0,lbg=lbg,ubg=ubg)
    with self.assertOutput(["Total number of inequality constraints...............:       10"],[]):
        ca.nlpsol("solver","ipopt",nlp,{"detect_simple_bounds": False})(x0=0,lbg=lbg,ubg=ubg)
    
    
    options = {"detect_simple_bounds": True,"ipopt.print_level":print_level,"print_time":False,"ipopt.tol":1e-12,"ipopt.fixed_variable_treatment":"relax_bounds","ipopt.bound_relax_factor":1e-12}
    solver = ca.nlpsol("solver","ipopt",nlp,options)
    options_nominal = {"detect_simple_bounds": False,"ipopt.print_level":print_level,"print_time":False,"ipopt.tol":1e-12,"ipopt.bound_relax_factor":1e-12}
    solver_nominal = ca.nlpsol("solver","ipopt",nlp,options_nominal)
    w_options = {"ipopt.warm_start_init_point":'yes',
                'ipopt.warm_start_bound_push':1e-16,
                "ipopt.bound_push":1e-16,
                "ipopt.warm_start_mult_bound_push": 1e-16,
                "ipopt.tol":1e-7}
    wsolver = ca.nlpsol("solver","ipopt",nlp,merge(options,w_options))
    wsolver_nominal = ca.nlpsol("solver","ipopt",nlp,merge(options_nominal,w_options))
    
    I = ca.densify(ca.DM.eye(5))
    
    for lbx,ubx in [(-ca.DM.inf(5),ca.DM.inf(5)),(-3*ca.DM.ones(5),3*ca.DM.ones(5)),(-4*ca.DM.ones(5),4*ca.DM.ones(5))]:
    
        for i in range(5):
            for sign in [-1,1]:
                x_target = I[:,i]*sign*20
                
                out = solver(p=x_target,lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
                out_nominal = solver_nominal(p=x_target,lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)

                if float(lbx[0])==-np.inf:
                    self.checkarray(out["lam_x"],ca.DM([0]*5))

                # Are the solutions equivalent
                self.checkarray(out["x"],out_nominal["x"],digits=8)
                self.checkarray(out["g"],out_nominal["g"],digits=8)
                
                # Can the compact output serve as zero-shot soluton of compact?
                wsolver(x0=out["x"],lam_x0=out["lam_x"],lam_g0=out["lam_g"],p=x_target,lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
                self.assertEqual(wsolver.stats()['iter_count'],0)
                
                # Can the compact output serve as zero-shot soluton of nominal?
                wsolver_nominal(x0=out["x"],lam_x0=out["lam_x"],lam_g0=out["lam_g"],p=x_target,lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
                self.assertEqual(wsolver_nominal.stats()['iter_count'],0)
                
                # Can the nominal output serve as zero-shot solution of compact?
                wsolver(x0=out_nominal["x"],lam_x0=out_nominal["lam_x"],lam_g0=out_nominal["lam_g"],p=x_target,lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
                self.assertEqual(wsolver.stats()['iter_count'],0)
 
  def test_derivative(self):
    x = ca.MX.sym("x",3)
    p = ca.MX.sym("p",3)
    nlp = {"x":x,"p":p,"f":ca.sumsqr(ca.sin(x)-p),"g":x}
    solver = ca.nlpsol("solver","sqpmethod",nlp,{"qpsol":"qrqp"})
    res = solver(p=p,x0=0,lbg=-inf,ubg=inf)
    f = ca.Function('f',[p],[res["x"]])
    f_ref = ca.Function('f_ref',[p],[ca.arcsin(p)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(0.1,0.2,0.3)])
    
    x = ca.MX.sym("x",3)
    p = ca.MX.sym("p",3)
    nlp = {"x":x,"p":p,"f":ca.sumsqr(ca.sin(x)-p)}
    solver = ca.nlpsol("solver","sqpmethod",nlp,{"qpsol":"qrqp"})
    res = solver(p=p,x0=0)
    f = ca.Function('f',[p],[res["x"]])
    f_ref = ca.Function('f_ref',[p],[ca.arcsin(p)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(0.1,0.2,0.3)])
    
  def test_hock_schittkowski(self):
    return
    import importlib
    file_names = []
    for i in range(1, 120):
        if i < 10:
            file_names.append("hs00" + str(i))
        elif i < 100:
            file_names.append("hs0" + str(i))
        else:
            file_names.append("hs" + str(i))




    solvers = []
    if ca.has_conic("qpoases"): solvers.append("qpoases")
    if ca.has_conic("qrqp"): solvers.append("qrqp")
    elastic_modes = [False, True]
    so_modes = [False, True]
    convex = "regularize"
    max_iter = 200
    
    import pandas as pd
    
    results = pd.DataFrame(columns=["problem","solver","elastic_mode","so_mode","outcome"])
    
    reference = pd.read_csv('hock_schittkowski/results.csv')
    


    for solver in solvers:
        for elastic_mode in elastic_modes:
            for so_mode in so_modes:
                outcomes_encountered = {"codegen": False}
                for file_name in file_names:
                    mode = (solver,elastic_mode,so_mode)
                    
                    if file_name == 'hs087':
                        continue
                    if file_name == 'hs082' or file_name == 'hs094' or file_name == 'hs115'\
                            or file_name == 'hs087':
                        continue

                    prob = importlib.import_module("hock_schittkowski."+file_name)
                    hock_schittkowsky_func = getattr(prob, file_name)
                    (x_opt,
                        f_opt, x, f, g, lbg, ubg, lbx, ubx, x0) = hock_schittkowsky_func()
                    # SOLVE NLP WITH CASADI-SQP
                    nlp = {'x': x, 'f': f, 'g': g}
                    e_o_f = False
                    with capture_stdout():
                        sqp_solver = ca.nlpsol('sqp_solver', 'sqpmethod', nlp, {"init_feasible":True, "elastic_mode": elastic_mode, 'qpsol': solver, "second_order_corrections": so_mode, 'convexify_strategy': convex, 'max_iter': max_iter, 'gamma_max': 10e40, 'qpsol_options.error_on_fail': e_o_f})
                    try:
                        with capture_stdout():
                            sol_sqp = sqp_solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
                            
                        outcome = sqp_solver.stats()["return_status"]

                    except Exception as ex:
                        outcome = "exception"
                    if outcome !="exception":
                        if outcome not in outcomes_encountered:
                            if outcome not in ["Search_Direction_Becomes_Too_Small","Non_Regular_Sensitivities"]: # Bugs
                                print("outcome:", outcome)
                                outcomes_encountered.add(outcome)
                                with capture_stdout():
                                    if solver=="qrqp": self.check_codegen(sqp_solver,dict(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg),std="c99",extra_options=["-Wno-maybe-uninitialized"])
                                
                    results = results.append({"problem":file_name,"solver":solver,"elastic_mode":elastic_mode,"so_mode": so_mode,"outcome": outcome},ignore_index=True)
                    
                    reference_outcome = list(reference[(reference['problem']==file_name) & (reference['solver']==solver) & (reference['elastic_mode']==elastic_mode) & (reference['so_mode']==so_mode)]["outcome"])[0]
                    self.assertEqual(outcome,reference_outcome)


    #results.to_csv('hock_schittkowski/results.csv',index=False)

  def test_exception_in_oraclefunction(self):
    x=ca.MX.sym("x")
    x_fail = x.attachAssert(x!=x,"Cuckoo")
    nlp={'x':x, 'f':(x-1)**2, 'g':x_fail}
        
    for Solver, solver_options, aux_options in solvers:
      if "madnlp" in Solver: continue
      
      print("test_exception_in_oraclefunction",Solver,solver_options)
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      with self.assertInAnyOutput("Cuckoo"):
        solver_out = solver(**solver_in)

  def test_unsolved_stats(self):
    x=ca.MX.sym("x")
    x_fail = x.attachAssert(x!=x,"Cuckoo")
    nlp={'x':x, 'f':(x-1)**2, 'g':x_fail}
        
    for Solver, solver_options, aux_options in solvers:
      solver = ca.nlpsol("mysolver", Solver, nlp, solver_options)
      with self.assertInException("No stats available"):
        solver.stats()

      
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipopt_detect_simple_bounds(self):
    x=ca.MX.sym("x")
    y=ca.MX.sym("y")
    z=ca.MX.sym("z")
    w=ca.MX.sym("w")
    
    def testcases():
    
        for d_scale in [1.0,1e5]:
            for c_scale in [1.0,1e5]:
    
                g = [x,(x+y)*d_scale,x**2+z**2,z,(z+y+w)*c_scale,z**2+y**2+w]
                lbg = [-5.0,-2*d_scale,-3.0,1.0,2*c_scale,3.0]
                ubg = [3.0,2*d_scale,3.0,1.0,2*c_scale,3.0]
                
                g = ca.vcat(g)
                lbg = ca.vcat(lbg)
                ubg = ca.vcat(ubg)
                
                lbx = [-inf,-inf,-inf,-inf]
                ubx = [inf,inf,inf,inf]
                lbx = ca.vcat(lbx)
                ubx = ca.vcat(ubx)
                
                yield (g,lbg,ubg,lbx,ubx,[0,3,5,1,4,2])
                
        
                g = [(x+y)*d_scale,x**2+z**2,(z+y+w)*c_scale,z**2+y**2+w]
                lbg = [-2*d_scale,-3.0,2*c_scale,3.0]
                ubg = [2*d_scale,3.0,2*c_scale,3.0]
                
                g = ca.vcat(g)
                lbg = ca.vcat(lbg)
                ubg = ca.vcat(ubg)
                
                lbx = [-5,-inf,1,-inf]
                ubx = [3,inf,1,inf]
                lbx = ca.vcat(lbx)
                ubx = ca.vcat(ubx)

                yield (g,lbg,ubg,lbx,ubx,[3,0,2,1])
    
    # 'make_parameter_nodual' was added in ipopt 3.14; an older ipopt rejects the
    # whole option set at construction ("Invalid options were detected by Ipopt"),
    # which has nothing to do with detect_simple_bounds.  Probe once.
    fixed_variable_treatments = ["make_constraint","make_parameter","make_parameter_nodual"]
    try:
        _p = ca.SX.sym("probe")
        ca.nlpsol("probe","ipopt",{"x":_p,"f":_p**2},
                  {"ipopt.fixed_variable_treatment":"make_parameter_nodual",
                   "ipopt.print_level":0,"print_time":False})
    except Exception:
        print("ipopt rejects fixed_variable_treatment=make_parameter_nodual "
              "(added in 3.14), skipping that variant")
        fixed_variable_treatments = fixed_variable_treatments[:-1]

    for (g,lbg,ubg,lbx,ubx,perm) in testcases():
        nlp = {"x":ca.vertcat(x,y,z,w),"g":g,"f":x**2+2*y**2+3*z**2+8*w**2}
        solver_ref = ca.nlpsol("solver","ipopt",nlp)
        ref = solver_ref(lbg=lbg,ubg=ubg,lbx=lbx,ubx=ubx)
        for g_perm in [range(g.numel()),perm]:
            nlp = {"x":ca.vertcat(x,y,z,w),"g":g[g_perm],"f":x**2+2*y**2+3*z**2+8*w**2}
            solver_ref2 = ca.nlpsol("solver","ipopt",nlp)
            ref2 = solver_ref2(lbg=lbg[g_perm],ubg=ubg[g_perm],lbx=lbx,ubx=ubx)
            for detect_simple_bounds in [True,False]:
                for fixed_variable_treatment in fixed_variable_treatments:
                    for start_with_resto in ["no","yes"]:
                        nlp = {"x":ca.vertcat(x,y,z,w),"g":g[g_perm],"f":x**2+2*y**2+3*z**2+8*w**2}
                        solver = ca.nlpsol("solver","ipopt",nlp,{"ipopt.fixed_variable_treatment":fixed_variable_treatment,"detect_simple_bounds": detect_simple_bounds,"ipopt.start_with_resto": start_with_resto})

                        res = solver(x0=ref["x"]*1.1,lbg=lbg[g_perm],ubg=ubg[g_perm],lbx=lbx,ubx=ubx)
                        print(g,detect_simple_bounds,g_perm,detect_simple_bounds,fixed_variable_treatment,start_with_resto)
                        for k in ["x","f"]:
                          self.checkarray(res[k],ref[k],digits=7)                        
                        for k in ["x","f","g"]+["lam_g"] if fixed_variable_treatment=="make_constraint" else []:
                          self.checkarray(res[k],ref2[k],digits=7)
                        self.assertTrue(solver.stats()["success"])
                        
                        cache = {"nlp_jac_g": solver.get_function("nlp_jac_g"), "nlp_hess_l": solver.get_function("nlp_hess_l")}
        
                        solver = ca.nlpsol("solver","ipopt",nlp,{"ipopt.fixed_variable_treatment":fixed_variable_treatment,"detect_simple_bounds": detect_simple_bounds,"ipopt.start_with_resto": start_with_resto,"cache":cache})

                        res = solver(x0=ref["x"]*1.1,lbg=lbg[g_perm],ubg=ubg[g_perm],lbx=lbx,ubx=ubx)
                        for k in ["f","x"]:
                          self.checkarray(res[k],ref[k],digits=7)
                        for k in ["x","f","g"]+["lam_g"] if fixed_variable_treatment=="make_constraint" else []:
                          self.checkarray(res[k],ref2[k],digits=7)
                          
                        # What if you pass a jacobian/hessian produced without regard for detect_simple_bounds?
                        cache = {"nlp_jac_g": solver_ref2.get_function("nlp_jac_g"), "nlp_hess_l": solver_ref2.get_function("nlp_hess_l")}       
        
                        solver = ca.nlpsol("solver","ipopt",nlp,{"ipopt.fixed_variable_treatment":fixed_variable_treatment,"detect_simple_bounds": detect_simple_bounds,"ipopt.start_with_resto": start_with_resto,"cache":cache})

                        res = solver(x0=ref["x"]*1.1,lbg=lbg[g_perm],ubg=ubg[g_perm],lbx=lbx,ubx=ubx)
                        for k in ["f","x"]:
                          self.checkarray(res[k],ref[k],digits=7)
                        for k in ["x","f","g"]+["lam_g"] if fixed_variable_treatment=="make_constraint" else []:
                          self.checkarray(res[k],ref2[k],digits=7)
                          
                        self.assertTrue(solver.stats()["success"])
        
  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipopt_custom_hess(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    w=ca.SX.sym("w")
    
    p = ca.SX(0,1)
    
    f = (1-x)**2+100*w**2
    g = y-x**2-w
    
    x = ca.vertcat(x,y,w)
    nlp={'x':x, 'f':f,'g': g}
    lam_f = ca.SX.sym("lam_f")
    lam_g = ca.SX.sym("lam_g",g.numel())
    
    ref_solver = ca.nlpsol("solver","ipopt",nlp)
    
    x0 = 0
    
    ref_sol = ref_solver(x0=0,lbg=0,ubg=0)
    
    
    nlp_hess_l_custom = ca.Function('nlp_hess_l',[x,p,lam_f,lam_g],[ca.triu(ca.DM.zeros(3,3))])
    options=  {}
    options["cache"] = {"nlp_hess_l":nlp_hess_l_custom}
    solver = ca.nlpsol("solver","ipopt",nlp,options)
    
    with self.assertInException("evaluation error"):
        self.checkfunction_light(ref_solver.get_function("nlp_hess_l"),solver.get_function("nlp_hess_l"),inputs=[ca.vertcat(0.11,0.3,0.7),0,1.2,3.7])
    

    lag = lam_f*f+ca.dot(lam_g,g)
    H = ca.jacobian(ca.gradient(lag,x),x,{"symmetric":True})

    
    nlp_hess_l_custom = ca.Function('nlp_hess_l',[x,p,lam_f,lam_g],[ca.triu(H)])
    
    ref_solver.get_function("nlp_hess_l").disp(True)
    nlp_hess_l_custom.disp(True)
    options=  {}
    options["cache"] = {"nlp_hess_l":nlp_hess_l_custom}
    solver = ca.nlpsol("solver","ipopt",nlp,options)
    
    sol = solver(x0=0,lbg=0,ubg=0)
    
    self.assertTrue(solver.stats()["success"])
    
    self.checkarray(sol["x"],ref_sol["x"],digits=6)
    
    self.checkfunction_light(ref_solver.get_function("nlp_hess_l"),solver.get_function("nlp_hess_l"),inputs=[ca.vertcat(0.11,0.3,0.7),0,1.2,3.7])
    

  @requires_nlpsol("ipopt")
  @memory_heavy()
  def test_ipopt_custom_jac(self):
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    w=ca.SX.sym("w")
    
    p = ca.SX(0,1)
    
    f = (1-x)**2+100*w**2
    g = y-x**2-w
    
    x = ca.vertcat(x,y,w)
    nlp={'x':x, 'f':f,'g': g}
    lam_f = ca.SX.sym("lam_f")
    lam_g = ca.SX.sym("lam_g",g.numel())
    
    ref_solver = ca.nlpsol("solver","ipopt",nlp)
    
    x0 = 0
    
    ref_sol = ref_solver(x0=0,lbg=0,ubg=0)
    
    options = {}
    solver = ca.nlpsol("solver","ipopt",nlp,options)

    nlp_jac_g_custom = ca.Function('nlp_jac_g',[x,p],[g,ca.jacobian(g,x)],["x","p"],["g","jac_g_x"])
    options["nlp_jac_g"] = {"nlp_jac_g":nlp_jac_g_custom}
    
    sol = solver(x0=0,lbg=0,ubg=0)
    self.assertTrue(solver.stats()["success"])
    self.checkfunction_light(ref_solver.get_function("nlp_jac_g"),solver.get_function("nlp_jac_g"),inputs=[ca.vertcat(0.11,0.3,0.7),0])
    
  @requires_nlpsol("ipopt")
  def test_option_propagation(self):
    x = ca.MX.sym('x')

    y = (x.printme(0)**2)
    options = {"common_options": {"helper_options": {"enable_fd":True,"enable_forward":False,"enable_reverse":False}}, "ipopt": {"resto.tol" :1e-8, "hessian_approximation":"limited-memory","max_iter":0}}
    print(options)
    solver = ca.nlpsol("solver","ipopt",{"x":x,"f":y},options)
    with capture_stdout() as result:
        solver(x0=1)

    self.assertTrue(len(result[1])==0)

    r = []
    for l in result[0].splitlines():
        if "|> 0 : " in l:
            r.append(float(l.split(":")[1]))
    print(r)
    spread = np.max(r)-np.min(r)
    print(spread)
    self.assertTrue(spread>0)
    self.assertTrue(spread<1e-4)
    

    """
    options = {"specific_options": {"nlp_hess_l": {"helper_options": {"enable_fd":True,"enable_forward":False,"enable_reverse":False}}}, "ipopt": {"resto.tol" :1e-8}}
    print(options)
    solver = nlpsol("solver","ipopt",{"x":x,"f":y},options)
    solver(x0=1)

    self.assertTrue(len(result[1])==0)

    r = []
    for l in result[0].splitlines():
        if "|> 0 : " in l:
            r.append(float(l.split(":")[1]))
    print(r)
    spread = np.max(r)-np.min(r)
    print(spread)
    self.assertTrue(spread>0)
    self.assertTrue(spread<1e-4)
    """

if __name__ == '__main__':
    unittest.main()
    print(solvers)
