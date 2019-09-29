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
import itertools

import os
#GlobalOptions.setCatchErrorsPython(False)

solvers= []

if has_nlpsol("worhp")  and not args.ignore_memory_heavy:
  solvers.append(("worhp",{"worhp": {"TolOpti":1e-9}},set()))
  #solvers.append(("worhp",{"TolOpti":1e-20,"TolFeas":1e-20,"UserHM": False}))
  pass

if has_nlpsol("ipopt"):
  solvers.append(("ipopt",{"print_time":False,"ipopt": {"tol": 1e-10, "derivative_test":"second-order","print_level":0}},set()))
  solvers.append(("ipopt",{"print_time":False,"ipopt": {"tol": 1e-10, "derivative_test":"first-order","hessian_approximation": "limited-memory","print_level":0}},set()))

if has_nlpsol("ipopt") and has_nlpsol("sqpmethod"):
  qpsol_options = {"nlpsol": "ipopt", "nlpsol_options": {"ipopt.tol": 1e-12,"ipopt.tiny_step_tol": 1e-20, "ipopt.fixed_variable_treatment":"make_constraint","ipopt.print_level":0,"print_time":False,"print_time":False} }
  solvers.append(("sqpmethod",{"qpsol": "nlpsol","qpsol_options": qpsol_options,"print_header":False,"print_iteration":False,"print_time":False},set()))
  solvers.append(("sqpmethod",{"qpsol": "nlpsol","qpsol_options": qpsol_options,"hessian_approximation": "limited-memory","tol_du":1e-10,"tol_pr":1e-10,"min_step_size":1e-14,"print_header":False,"print_iteration":False,"print_time":False},set()))

if has_conic("qrqp") and has_nlpsol("sqpmethod"):
  qpsol_options = {"print_iter":False,"print_header":False,"error_on_fail" : False}
  solvers.append(("sqpmethod",{"qpsol": "qrqp","qpsol_options": qpsol_options,"print_header":False,"print_iteration":False,"print_time":False},{"codegen"}))
  solvers.append(("sqpmethod",{"qpsol": "qrqp","max_iter_ls":0,"qpsol_options": qpsol_options,"print_header":False,"print_iteration":False,"print_time":False},{"codegen"}))
  solvers.append(("sqpmethod",{"qpsol": "qrqp","convexify_strategy":"regularize","max_iter":500,"qpsol_options": qpsol_options,"print_header":False,"print_iteration":True,"print_time":False,"tol_du":1e-8,"min_step_size":1e-12},{"codegen"}))

if has_nlpsol("blocksqp"):
  try:
    load_linsol("ma27")
    solvers.append(("blocksqp",{},set()))
  except:
    pass

if has_nlpsol("bonmin"):
  solvers.append(("bonmin",{},{"discrete"}))

if "SKIP_KNITRO_TESTS" not in os.environ and has_nlpsol("knitro"):
  solvers.append(("knitro",{"knitro":{"feastol":1e-8,"opttol":1e-8}},set()))

print(solvers)

class NLPtests(casadiTestCase):

  @memory_heavy()
  def test_nonregular_point(self):
    x=SX.sym("x")

    nlp={'x':x,'f':(x+1)**2, 'g': sqrt(x)}

    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if Solver not in ["ipopt","snopt","blocksqp","bonmin"]:
        self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_NAN")
      self.assertFalse(solver.stats()["success"])

    nlp={'x':x,'f':x**2, 'g': sqrt(x)}
    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if Solver not in ["ipopt","snopt","bonmin"]:
        self.assertTrue(solver.stats()["unified_return_status"]=="SOLVER_RET_NAN")
      self.assertFalse(solver.stats()["success"])

  def test_iteration_interrupt(self):
   for Solver, solver_options, features in solvers:
      if Solver not in ["ipopt","sqpmethod"]: continue

      opti = Opti()

      x = opti.variable()

      def interrupt(i):
        raise KeyboardInterrupt()

      opti.minimize((x-1)**4)


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
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1.4-x)**2+100*(y-x**2)**2}

    for Solver, solver_options, features in solvers:
      if "discrete" not in features: continue
      self.message(str(Solver))
      solver_options = dict(solver_options)
      solver_options["discrete"] = [1,0]
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0.16,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))

      self.check_serialize(solver,solver_in)

    for Solver, solver_options, features in solvers:
      if "discrete" not in features: continue
      self.message(str(Solver))
      solver_options = dict(solver_options)
      solver_options["discrete"] = [0,1]
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],2e-4,3,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],sqrt(2),3,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],2,3,str(Solver))

      self.check_serialize(solver,solver_in)
  def test_nan(self):
    x=SX.sym("x")
    nlp={'x':x, 'f':-x,'g':x}
    for Solver, nlp_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, nlp_options)

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

  def test_wrongdims(self):
    x=SX.sym("x",2)
    nlp={'x':x, 'f':-x[0],'g':diag(x)}

    for Solver, solver_options, features in solvers:
      with self.assertInException("dense vector"):
        solver = nlpsol("mysolver", Solver, nlp, solver_options)
    nlp={'x':x, 'f':-x[0],'g':mtimes(x,x.T)}

    for Solver, solver_options, features in solvers:
      with self.assertInException("dense vector"):
        solver = nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'f':SX(1,1),'g':x}

    for Solver, solver_options, features in solvers:
      with self.assertInException("dense"):
        solver = nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'f':SX.zeros(0,0),'g':x}

    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)

    nlp={'x':x, 'g':x}
    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)

    x = vec(diag(SX.sym("x",2)))
    nlp={'x':x, 'f':mtimes(x.T,x),'g':x[0]}
    for Solver, solver_options, features in solvers:
      with self.assertInException("dense vector"):
        solver = nlpsol("mysolver", Solver, nlp, solver_options)


  def test_initialcond(self):
    x=SX.sym("x")
    nlp={'x':x, 'f':-cos(x),'g':x}

    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[-20]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      with self.assertRaises(Exception):
        solver_out = solver(**solver_in)

    for Solver, solver_options, features in solvers:
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[-20]
      with self.assertRaises(Exception):
        solver_out = solver(**solver_in)

  def test_IPOPT(self):
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, features in solvers:
      self.message("trivial " + str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPT_par(self):
    x=SX.sym("x")
    p=SX.sym("p")
    nlp={'x':x, 'p':p, 'f':(x-p)**2, 'g':x}

    for Solver, solver_options, features in solvers:
      self.message("trivial " + str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTinf(self):
    self.message("trivial IPOPT, infinity bounds")
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':x}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-inf]
      solver_in["ubx"]=[inf]
      solver_in["lbg"]=[-inf]
      solver_in["ubg"]=[inf]

      if Solver in ("worhp","knitro"):
        with self.assertRaises(Exception):
          solver_out = solver(**solver_in)
        return




      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver) + str(solver_out["x"][0]-1))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      if Solver in ("worhp","knitro"):
        with self.assertRaises(Exception):
          solver_out = solver(**solver_in)
        return

      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver) + str(solver_out["x"][0]-1))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,9,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,9,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")



  def test_IPOPTrb(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2}

    for Solver, solver_options, features in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTrb2(self):
    self.message("rosenbrock, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertTrue(solver.stats()["success"])

      digits = 6

      self.assertAlmostEqual(solver_out["f"][0],0,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0,5,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTrbf(self):
    self.message("rosenbrock fixed, limited-memory hessian approx")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_warmstart(self):

    x=SX.sym("x")
    y=SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    for Solver, solver_options, features in solvers:
      self.message(Solver)
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

      self.message(":warmstart")
      if "ipopt" in str(Solver):
        oldsolver=solver
        options = dict(solver_options)
        options["ipopt.warm_start_init_point"]="yes"
        options["ipopt.warm_start_bound_push"]=1e-6
        options["ipopt.warm_start_slack_bound_push"]=1e-6
        options["ipopt.warm_start_mult_bound_push"]=1e-6
        options["ipopt.mu_init"]=1e-6
        solver = nlpsol("mysolver", Solver, nlp, options)

        solver_in["lbx"]=[-10]*2
        solver_in["ubx"]=[10]*2
        solver_in["lbg"]=[0]
        solver_in["ubg"]=[1]
        solver_in["x0"]=oldsolver_out["x"]
        solver_in["lam_g0"]=oldsolver_out["lam_g"]
        if "bonmin" not in str(Solver): solver_in["lam_x0"] =oldsolver_out["lam_x"]


        solver_out = solver(**solver_in)

  def test_IPOPTrhb2_gen(self):
    self.message("rosenbrock, exact hessian generated, constrained")
    x=SX.sym("x")
    y=SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    sigma=SX.sym("sigma")
    lambd=SX.sym("lambd")

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))
      self.check_serialize(solver, solver_in)

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")


  def test_jacG_empty(self):
    x=SX.sym("x")
    y=SX.sym("y")

    obj = (1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'f':obj, 'g':1}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      if "worhp"==Solver:
        continue
      if "sqpmethod"==Solver:
        continue
      if "snopt"==Solver:
        continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[2]
      solver_out = solver(**solver_in)

      digits = 5

      self.checkarray(solver_out["f"],DM([0]),str(Solver),digits=digits)
      self.checkarray(solver_out["x"],DM([1,1]),str(Solver),digits=digits)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0,0]),str(Solver),digits=digits)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],DM([0]),str(Solver),digits=digits)

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTrhb2_gen_par(self):
    self.message("rosenbrock, exact hessian generated, constrained, parametric")
    x=SX.sym("x")
    y=SX.sym("y")
    p=SX.sym("p")

    obj = (p-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'p':p, 'f':obj, 'g':x**2+y**2}

    c_r = 4.56748075136258e-02;
    x_r = [7.86415156987791e-01,6.17698316967954e-01]

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[0.5,0.5]
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["lbg"]=[0]
      solver_in["ubg"]=[1]
      solver_in["p"]=[1]


      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
      solver_out = solver(**solver_in)

      digits = 5

      self.assertAlmostEqual(solver_out["f"][0],c_r,digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],x_r[0],digits,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],x_r[1],digits,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,5 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,5 if Solver=="snopt" else 8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][0],0.12149655447670,6,str(Solver))


  def test_IPOPTrhb_gen(self):
    self.message("rosenbrock, exact hessian generated")
    x=SX.sym("x")
    y=SX.sym("y")

    obj=(1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'f':obj}

    sigma=SX.sym("sigma")

    for Solver, solver_options, features in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,7,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,8,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,8,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTrhb_gen_xnonfree(self):
    self.message("rosenbrock, exact hessian generated, non-free x")
    x=SX.sym("x")
    y=SX.sym("y")

    obj=(1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'f':obj}

    sigma=SX.sym("sigma")

    for Solver, solver_options, features in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[1,-10]
      solver_in["ubx"]=[1,10]



      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTrhb_gen_par(self):
    self.message("rosenbrock, exact hessian generated, parametric")
    x=SX.sym("x")
    y=SX.sym("y")

    p=SX.sym("p")
    obj=(p-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(*[x,y]), 'p':p, 'f':obj}

    sigma=SX.sym("sigma")

    for Solver, solver_options, features in solvers:
      if "sqpmethod"==Solver and "regularize" in str(solver_options): continue
      if "snopt"==Solver: continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]*2
      solver_in["ubx"]=[10]*2
      solver_in["p"]=1
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1,7,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
  @memory_heavy()
  def test_IPOPTnorm(self):
    self.message("IPOPT min ||x||^2_2")
    def norm_2(mx):
      return c.dot(mx,mx)
    N=10
    x=MX.sym("x",N)
    x0=linspace(0,1,N)
    X0=MX(x0)
    nlp={'x':x, 'f':norm_2(x-X0), 'g':2*x}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0]*10),8,str(Solver),digits=8)
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_g"][1],0,8,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
  def test_IPOPTnoc(self):
    self.message("trivial IPOPT, no constraints")
    """ There is an assertion error thrown, but still it works"""
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2}
    for Solver, solver_options, features in solvers:
      if "snopt"==Solver: continue
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      # ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,7,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTmx(self):
    self.message("trivial IPOPT, using MX")
    x=MX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':2*x}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      # ({"tol":1e-8,"max_iter":103, "MaxIter": 103,"print_level":0,"derivative_test":"first-order"}).iteritems():

      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10]
      solver_in["ubg"]=[10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,9,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTc(self):
    self.message("trivial, overconstrained")
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':vertcat(*[x,x,x])}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10, -10, -10]
      solver_in["ubg"]=[10, 10, 10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,9,str(Solver) )
      self.assertAlmostEqual(solver_out["x"][0],1,5,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTc2(self):
    self.message("trivial2, overconstrained")
    x=SX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':vertcat(*[x,x,x+x])}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10, -10, -10]
      solver_in["ubg"]=[10, 10, 10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,10,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,8,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
  def test_IPOPTcmx(self):
    self.message("trivial , overconstrained, using MX")
    x=MX.sym("x")
    nlp={'x':x, 'f':(x-1)**2, 'g':vertcat(*[2*x,3*x,4*x])}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10]
      solver_in["ubx"]=[10]
      solver_in["lbg"]=[-10,-10,-10]
      solver_in["ubg"]=[10,10,10]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["f"][0],0,9,str(Solver))
      self.assertAlmostEqual(solver_out["x"][0],1,8,str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTdeg(self):
    self.message("degenerate optimization IPOPT")
    x=SX.sym("x")
    y=SX.sym("y")
    nlp={'x':vertcat(*[x,y]), 'f':0, 'g':vertcat(*[x-y,x])}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10, -10]
      solver_in["ubx"]=[10, 10]
      solver_in["lbg"]=[0, 3]
      solver_in["ubg"]=[0, 3]
      solver_out = solver(**solver_in)
      self.assertAlmostEqual(solver_out["x"][0],solver_out["x"][1],4 if "sqic" in str(solver_options) else 10,"IPOPT")

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_IPOPTdegc(self):
    self.message("degenerate optimization IPOPT, overconstrained")
    x=SX.sym("x")
    y=SX.sym("y")
    nlp={'x':vertcat(*[x,y]), 'f':0, 'g':vertcat(*[x-y,x,x+y])}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=[-10, -10]
      solver_in["ubx"]=[10, 10]
      solver_in["lbg"]=[0, 3 , -10]
      solver_in["ubg"]=[0, 3, 10]
      solver_out = solver(**solver_in)
      # todo: catch error when set([0, 3 , 5]) two times
      self.assertAlmostEqual(solver_out["x"][0],solver_out["x"][1],4 if "sqic" in str(solver_options) else 10,"IPOPT")

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_XfreeChange(self):
    self.message("Change in X settings")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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
      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_activeLBX(self):
    self.message("active LBX")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(Solver)
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
  def test_activeLBG(self):
    self.message("active LBG")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_activeUBG(self):
    self.message("active UBG")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_activeUBX(self):
    self.message("active UBX")
    x=SX.sym("x")
    y=SX.sym("y")

    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+100*(y-x**2)**2, 'g':x+y}
    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
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

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")
  @memory_heavy()
  def test_QP(self):
    self.message("QP")

    N = 50

    x = SX.sym("x",N)
    x0 = DM(list(range(N)))
    H = diag(list(range(1,N+1)))
    obj = 0.5*mtimes([(x-x0).T,H,(x-x0)])

    nlp = {'x':x, 'f':obj}
    for Solver, solver_options, features in solvers:
      if "snopt"==Solver: continue
      self.message(str(Solver))
      if Solver=="sqpmethod" and "limited-memory" in str(solver_options): continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}
      solver_in["lbx"]=-1000
      solver_in["ubx"]=1000
      solver_out = solver(**solver_in)
      self.checkarray(solver_out["x"],x0,str(Solver),digits=2)
      self.assertAlmostEqual(solver_out["f"][0],0,3,str(Solver))
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM.zeros(N,1),str(Solver),digits=4)

  def test_QP2(self):
    H = DM([[1,-1],[-1,2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])

    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0.5,0])
    UBX = DM([0.5,inf])

    x=SX.sym("x",2)
    nlp={'x':x, 'f':0.5*mtimes([x.T,H,x])+mtimes(G.T,x), 'g':mtimes(A,x)}

    for Solver, solver_options, features in solvers:
      self.message(Solver)
      options = dict(solver_options)
      if "ipopt" in str(Solver):
        options["ipopt.fixed_variable_treatment"] = "make_constraint"
      solver = nlpsol("mysolver", Solver, nlp, options)
      #{"tol":1e-8,"TolOpti":1e-25,"hessian_approximation":"limited-memory","max_iter":100,"MaxIter": 100,"print_level":0, "fixed_variable_treatment": "make_constraint"}).iteritems():
      solver_in = {}

      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA
      if 'sqic' in str(solver_options):
        continue

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1.25,6,str(Solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],4.75,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],DM([0,2,0]),str(Solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-7.4375,6,str(Solver))
      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

      solver = nlpsol("mysolver", Solver, nlp, options)
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],0.5,6,str(Solver))
      self.assertAlmostEqual(solver_out["x"][1],1.25,6,str(Solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],4.75,6,str(Solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(Solver))

      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_g"],DM([0,2,0]),str(Solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-7.4375,6,str(Solver))
      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_QP2_unconvex(self):
    H = DM([[1,-1],[-1,-2]])
    G = DM([-2,-6])
    A =  DM([[1, 1],[-1, 2],[2, 1]])

    LBA = DM([-inf]*3)
    UBA = DM([2, 2, 3])

    LBX = DM([0]*2)
    UBX = DM([inf]*2)

    x=SX.sym("x",2)
    nlp={'x':x, 'f':0.5*mtimes([x.T,H,x])+mtimes(G.T,x), 'g':mtimes(A,x)}

    for Solver, solver_options, features in solvers:
      self.message(Solver)
      options = dict(solver_options)
      if "ipopt" in str(Solver):
        options["ipopt.fixed_variable_treatment"] = "make_constraint"
      solver = nlpsol("mysolver", Solver, nlp, options)
      solver_in = {}
      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      if "qrqp" in str(solver_options):
        solver_in["x0"]= DM([1,1])
        solver_in["lam_g0"]= DM([1,1,0])
      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,6,str(solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(solver))

      if "bonmin" not in str(Solver) and "worhp" not in str(Solver): self.checkarray(solver_out["lam_g"],DM([4+8.0/9,20.0/9,0]),str(solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-10-16.0/9,6,str(solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

      solver = nlpsol("mysolver", Solver, nlp, options)

      solver_in["lbx"]=LBX
      solver_in["ubx"]=UBX
      solver_in["lbg"]=LBA
      solver_in["ubg"]=UBA

      solver_out = solver(**solver_in)

      self.assertAlmostEqual(solver_out["x"][0],2.0/3,6,str(solver))
      self.assertAlmostEqual(solver_out["x"][1],4.0/3,6,str(solver))

      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][0],0,6,str(solver))
      if "bonmin" not in str(Solver): self.assertAlmostEqual(solver_out["lam_x"][1],0,6,str(solver))

      if "bonmin" not in str(Solver) and "worhp" not in str(Solver): self.checkarray(solver_out["lam_g"],DM([4+8.0/9,20.0/9,0]),str(solver),digits=6)

      self.assertAlmostEqual(solver_out["f"][0],-10-16.0/9,6,str(solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_bug(self):
    x = MX.sym("x", 3)
    y = MX.sym("y", 2)
    f = Function("f", [x, y], [1.])

    aa = MX.sym("aa", 5)
    a = aa[:3]
    b = aa[3:]
    f_call = f(a, b)
    nlp = {'x':aa, 'f':f_call}
    for Solver, solver_options, features in solvers:
      if "worhp" in Solver: continue
      if "snopt"==Solver: continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

  def test_missing_symbols(self):
    x = MX.sym("x")
    p = MX.sym("p")

    for Solver, solver_options, features in solvers:
      with self.assertInException("[p] are free"):
        solver = nlpsol("solver",Solver,{"x":x,"f":(x-p)**2}, solver_options)

  @requires_nlpsol("ipopt")
  def test_no_success(self):

    x=SX.sym("x")
    y=SX.sym("y")

    f = (1-x)**2+100*(y-x**2)**2
    for Solver, solver_options, features in solvers:
      solver = nlpsol("solver","ipopt",{'x':vertcat(x,y), 'f':f,'g':vertcat(x+1,x-2)})
      solver(x0=0,lbg=0,ubg=0)
      self.assertFalse(solver.stats()["success"])

      solver = nlpsol("solver","ipopt",{'x':vertcat(x,y), 'f':f,'g':vertcat(x+1,x-2)},{"error_on_fail":True})
      with self.assertInException("process"):
        solver(x0=0,lbg=0,ubg=0)


  @requires_nlpsol("ipopt")
  def test_iteration_Callback(self):

    x=SX.sym("x")
    y=SX.sym("y")

    f = (1-x)**2+100*(y-x**2)**2
    nlp={'x':vertcat(x,y), 'f':f,'g':x+y}
    fcn = Function('f', [x, y], [f])

    class MyCallback(Callback):
      def __init__(self,nx, ng, np):
        Callback.__init__(self)
        self.foo = []

        self.nx = nx
        self.ng = ng
        self.np = np
        self.construct("mycallback", {})

      def get_n_in(self): return nlpsol_n_out()
      def get_n_out(self): return 1


      def get_sparsity_in(self, i):
        n = nlpsol_out(i)
        if n=='f':
          return Sparsity. scalar()
        elif n in ('x', 'lam_x'):
          return Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
          return Sparsity.dense(self.ng)
        else:
          return Sparsity(0,0)
      def eval(self, arg):
        self.foo.append(arg)
        return [0]

    mycallback = MyCallback(2,1,0)
    opts = {}
    opts['iteration_callback'] = mycallback
    opts['ipopt.tol'] = 1e-8
    opts['ipopt.max_iter'] = 50
    solver = nlpsol('solver', 'ipopt', nlp, opts)
    sol = solver(lbx=-10, ubx=10, lbg=-10, ubg=10)
    self.assertEqual(len(mycallback.foo),solver.stats()["iter_count"]+1)

    class MyCallback(Callback):
      def __init__(self,nx, ng, np):
        Callback.__init__(self)
        self.foo = []

        self.nx = nx
        self.ng = ng
        self.np = np
        self.construct("mycallback", {})

      def get_n_in(self): return nlpsol_n_out()
      def get_n_out(self): return 1


      def get_sparsity_in(self, i):
        n = nlpsol_out(i)
        if n=='f':
          return Sparsity. scalar()
        elif n in ('x', 'lam_x'):
          return Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
          return Sparsity.dense(self.nx)
        else:
          return Sparsity(0,0)
      def eval(self, arg):
        self.foo.append(arg)
        return [0]

    mycallback = MyCallback(2,1,0)
    opts = {}
    opts['iteration_callback'] = mycallback
    opts['ipopt.tol'] = 1e-8
    opts['ipopt.max_iter'] = 50

    try:
      solver = nlpsol('solver', 'ipopt', nlp, opts)
    except Exception as e:
      self.assertTrue("Callback function input size mismatch" in str(e))
  def test_pathological(self):
    x=SX.sym("x")
    y=SX.sym("y")
    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+y**2}

    for Solver, solver_options, features in solvers:
      if "snopt"==Solver: continue
      self.message(str(Solver))
      if "worhp"==Solver or "stabilizedsqp"==Solver : continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,-1]
      solver_in["ubx"]=[10,2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],DM([0]),digits=7)
      self.checkarray(solver_out["x"],DM([1,0]),digits=7,failmessage=str(Solver))
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0,-0]),digits=7,failmessage=str(Solver))

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_pathological2(self):
    x=SX.sym("x")
    y=SX.sym("y")
    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2+y}

    for Solver, solver_options, features in solvers:
      if "snopt"==Solver: continue
      self.message(Solver)
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,0]
      solver_in["ubx"]=[10,2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],DM([0]),digits=7)
      self.checkarray(solver_out["x"],DM([1,0]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0,-1]),digits=7)

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_pathological3(self):
    x=SX.sym("x")
    y=SX.sym("y")
    nlp={'x':vertcat(*[x,y]), 'f':(1-x)**2, 'g':x+y}

    for Solver, solver_options, features in solvers:
      self.message(str(Solver))
      if "worhp"==Solver: continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[1,1]
      solver_in["lbx"]=[-10,0]
      solver_in["ubx"]=[10,2]
      solver_in["lbg"]=[2]
      solver_in["ubg"]=[2]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],DM([0]),digits=7)
      self.checkarray(solver_out["x"],DM([1,1]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0,0]),digits=7)

      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_pathological4(self):
    x=SX.sym("x")
    nlp={'x':x, 'f':x*x}

    for Solver, solver_options, features in solvers:
      if "snopt"==Solver: continue
      self.message(Solver)
      if "worhp"==Solver: continue
      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      solver_in = {}

      solver_in["x0"]=[0]
      solver_in["lbx"]=[0]
      solver_in["ubx"]=[0]

      solver_out = solver(**solver_in)

      self.checkarray(solver_out["f"],DM([0]),digits=7)
      self.checkarray(solver_out["x"],DM([0]),digits=7)
      if "bonmin" not in str(Solver): self.checkarray(solver_out["lam_x"],DM([0]),digits=7)
      if "codegen" in features:
        self.check_codegen(solver,solver_in,std="c99")

  def test_nlp_sensitivity(self):

    x = MX.sym("x")
    p = MX.sym("p")

    nlp = {"x":x,"p":p,"f":(sin(x)-p**2)**2,"g":x}

    for Solver, solver_options, features in solvers:
      if "ipopt" in str(solver_options): continue

      solver = nlpsol("mysolver", Solver, nlp, solver_options)

      z = solver(p=p,x0=x,lbg=0)["x"]

      z2 = asin(p**2)

      f = Function('f',[x,p],[z,jacobian(z,p)])
      f2 = Function('f',[x,p],[z2,jacobian(z2,p)])

      self.checkfunction_light(f,f2,[0,0.5],digits=6)

  @requires_conic("qrqp")
  def test_regularize_sqpmethod(self):

    # Test problem that is indefinite in direction of the constraint Jacobian
    x = MX.sym("x",2)
    f = 0.5*bilin(DM([[1,0],[0,-2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    solver = nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True}})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats = solver.stats()
    self.assertTrue(stats["iter_count"]==1)
    self.assertTrue("H:\n[[1, 0], \n [0, -2]]" in result[0])
    self.checkarray(res["x"],DM([0,2]),digits=6)

    solver = nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize","convexify_margin":0})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==2)
    self.assertTrue("H:\n[[3, 0], \n [0, 0]]" in result[0])

    solver = nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize","convexify_margin":1e-4})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==2)
    self.assertTrue("H:\n[[3.0001, 0], \n [0, 0.0001]]" in result[0])

    x = MX.sym("x",2)
    f = 0.5*bilin(DM([[1,0],[0,2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    solver = nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True}})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats = solver.stats()
    self.assertTrue(stats["iter_count"]==1)
    print(result[0])
    self.assertTrue("H:\n[[1, 0], \n [0, 2]]" in result[0])
    self.checkarray(res["x"],DM([0,2]),digits=6)

    solver = nlpsol("mysolver", "sqpmethod", nlp, {"qpsol":"qrqp","qpsol_options": {"print_problem":True},"convexify_strategy":"regularize"})
    with capture_stdout() as result:
      res = solver(lbg=2,ubg=2)
    stats_reg = solver.stats()
    self.checkarray(res["x"],DM([0,2]),digits=6)
    self.assertTrue(stats_reg["iter_count"]==1)
    self.assertTrue("H:\n[[1, 0], \n [0, 2]]" in result[0])

  def test_indefinite(self):

    # Test problem that is indefinite in direction of the constraint Jacobian
    x = MX.sym("x",2)
    f = 0.5*bilin(DM([[1,0],[0,-2]]),x,x)

    nlp = {"x":x,"f":f,"g":x[1]}

    for Solver, solver_options, features in solvers:
      solver_in = {"lbg": 2, "ubg": 2}

      solver = nlpsol("mysolver", Solver, nlp, solver_options)
      out = solver(**solver_in)
      self.checkarray(out["x"],DM([0,2]),digits=6)
      if "bonmin" not in str(Solver): self.checkarray(out["lam_g"],DM([4]),digits=6)

      if "codegen" in features:
        solver.generate('f.c',{"main":True})
        solver.generate_in("in.dat",solver.convert_in(solver_in))
        print(solver_in)
        self.check_codegen(solver,solver_in,std="c99")

  @requires_nlpsol("sqpmethod")
  def test_gauss_newton_sqpmethod(self):
    x = SX.sym("x",3)

    F = sin(x) - vertcat(1,2,3)*vertcat(x[2],0,0)
    J = jacobian(F,x)
    f = 0.5*dot(F,F)
    p = SX.sym("x",0,1)
    lam_f = SX.sym("x")
    lam_g = SX.sym("x",0,1)
    GN = Function('GN',[x,p,lam_f,lam_g],[lam_f*triu(mtimes(J.T,J))])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    with self.assertInException("Hessian must be symmetric"):
      solver = nlpsol("solver","sqpmethod",nlp,options)

    # A 2-norrm problem ...
    F = sin(x-vertcat(0.22,0.72,0.2)) - vertcat(0.1,0.5,0.99)
    J = jacobian(F,x)
    f = 0.5*dot(F,F)
    H = Function("H",[x],[hessian(f,x)[0]])
    H = H(0)
    # with an indefinite Hessian at x0
    self.assertTrue(np.any(np.linalg.eig(H)[0]<0))

    # Solve with Gauss-Newton -> 6 iterations
    GN = Function('GN',[x,p,lam_f,lam_g],[lam_f*mtimes(J.T,J)])
    options = {"convexify_strategy":"regularize","qpsol":"qrqp","hess_lag": GN}
    nlp = {"x":x,"f":f}
    solver = nlpsol("solver","sqpmethod",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==6)

    # Solve with exact Hessian + regularization -> 9 iterations
    options = {"convexify_strategy":"regularize","qpsol":"qrqp"}
    nlp = {"x":x,"f":f}
    solver = nlpsol("solver","sqpmethod",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==9)

  @requires_nlpsol("ipopt")
  def test_gauss_newton_ipopt(self):
    x = SX.sym("x",3)

    F = sin(x) - vertcat(1,2,3)*vertcat(x[2],0,0)
    J = jacobian(F,x)
    f = 0.5*dot(F,F)
    p = SX.sym("x",0,1)
    lam_f = SX.sym("x")
    lam_g = SX.sym("x",0,1)
    GN = Function('GN',[x,p,lam_f,lam_g],[lam_f*mtimes(J.T,J)])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    with self.assertInException("Hessian must be upper triangular"):
      solver = nlpsol("solver","ipopt",nlp,options)

    # A 2-norrm problem ...
    F = sin(x-vertcat(0.22,0.72,0.2)) - vertcat(0.1,0.5,0.99)
    J = jacobian(F,x)
    f = 0.5*dot(F,F)
    H = Function("H",[x],[hessian(f,x)[0]])
    H = H(0)
    # with an indefinite Hessian at x0
    self.assertTrue(np.any(np.linalg.eig(H)[0]<0))

    # Solve with Gauss-Newton -> 6 iterations
    GN = Function('GN',[x,p,lam_f,lam_g],[lam_f*triu(mtimes(J.T,J))])
    options = {"hess_lag": GN}
    nlp = {"x":x,"f":f}
    solver = nlpsol("solver","ipopt",nlp,options)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==6)

    # Solve with exact Hessian + regularization -> 9 iterations
    nlp = {"x":x,"f":f}
    solver = nlpsol("solver","ipopt",nlp)
    res = solver()
    stats_reg = solver.stats()
    self.assertTrue(stats_reg["iter_count"]==9)

  def test_cvx_sqpmethod(self):

    eps = 1e-2

    H = diagcat(DM([[2,3],[3,0]]),DM([[7,1,0],[1,0,6],[0,6,3]]))
    p = [0,3,2,1,4]
    Hp = H[p,p]

    eig = np.linalg.eig

    def reflect(A,eps):
      [d,V] = np.linalg.eig(A)  
      d = abs(d)
      D = diag(d)
      return np.dot(V,np.dot(D,V.T))

    def clip(A,eps):
      [d,V] = np.linalg.eig(A)  
      d[d<eps] = eps
      D = diag(d)
      return np.dot(V,np.dot(D,V.T))


    for H, Hcvx, opts, includes_init, excludes_init in [
      (sparsify(DM([[0, 0],[0, 2]])), DM([[eps, 0],[0, 2+eps]]), {"convexify_strategy": "regularize","convexify_margin":eps},[],[]),
      (diagcat(DM([[2, 0,0],[0, -1,0],[0,0,8]]),DM([[5, 3],[3, 6]])), diagcat(DM([[2, 0,0],[0, 1,0],[0,0,8]]),DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 2 blocks with maximum size 3"],[]),
      (diagcat(DM([[2, 0,0],[0, -1,0],[0,0,8]]),DM([[5, 3],[3, 6]])), diagcat(DM([[2, 0,0],[0, eps,0],[0,0,8]]),DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      (DM([[2, 0,0],[0, eps/2,0],[0,0,8]]), DM([[2, 0,0],[0, eps,0],[0,0,8]]), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (diagcat(DM([[2, 0,0],[0, -1,0],[0,0,8]]),DM([[5, 8],[8, 6]])), diagcat(DM([[2, 0,0],[0, 1,0],[0,0,8]]),reflect(DM([[5, 8],[8, 6]]),eps)), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (diagcat(DM([[2, 0,0],[0, -1,0],[0,0,8]]),DM([[5, 8],[8, 6]])), diagcat(DM([[2, 0,0],[0, eps,0],[0,0,8]]),clip(DM([[5, 8],[8, 6]]),eps)), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      (DM([[2, 0,0],[0, eps/2,0],[0,0,8]]), DM([[2, 0,0],[0, eps,0],[0,0,8]]), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (diagcat(sparsify(DM([[2, 0,0],[0, -1,0],[0,0,8]])),DM([[5, 3],[3, 6]])), diagcat(DM([[2, 0,0],[0, 1,0],[0,0,8]]),DM([[5, 3],[3, 6]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 4 blocks with maximum size 2"],[]),
      (diagcat(sparsify(DM([[2, 0,0],[0, -1,0],[0,0,8]])),DM([[5, 3],[3, 5]])), diagcat(DM([[2, 0,0],[0, 1,0],[0,0,8]]),DM([[5, 3],[3, 5]])), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},["Identified 4 blocks with maximum size 2"],[]),
      (Hp, reflect(Hp,eps), {"convexify_strategy": "eigen-reflect","convexify_margin":eps},[],[]),
      (Hp, clip(Hp,eps), {"convexify_strategy": "eigen-clip","convexify_margin":eps},[],[]),
      ]:

      n = H.shape[0]
      x = MX.sym("H",n)
      nlp = {"x": x, "f": 0.5*bilin(H,x,x)}
      options = {"max_iter":1,"qpsol":"qrqp","verbose":True}
      options.update(opts)
      with self.assertOutput(includes_init,excludes_init): 
        solver = nlpsol("solver","sqpmethod",nlp,options)

      x0 = DM.ones(n,1)

      res = solver(x0=x0)

      self.checkarray(x0-np.linalg.solve(Hcvx,mtimes(H,x0)),res["x"])

      self.check_serialize(solver,{"x0":x0})

      self.check_codegen(solver,{"x0":x0},std="c99")


  def test_simple_bounds_detect(self):

    x = SX.sym("x",5)
    p = SX.sym("p",1)


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

    lbg = vcat(lbg)
    ubg = vcat(ubg)
    g = vcat(g)

    [gi,lbx,ubx,lam_f,lam_b]=detect_simple_bounds(x,p,g,lbg,ubg)

    self.checkarray(DM(lbx).T,[-4,9,7,-inf,0])
    self.checkarray(DM(ubx).T,[2,9,7,inf,2])

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
      a = DM.zeros(12,1)
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
      a = DM.zeros(12,1)
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
      a = DM.zeros(5,1)
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
      a = DM.zeros(5,1)
      a[i] = 2
      b = round_trip_b(([3,7],a))
      c = round_trip_b(round_trip_b(([3,7],a)))
      self.checkarray(b[0].T,[3,7])
      self.checkarray(b[1],G[i,:])
      self.checkarray(b[1],c[1])

if __name__ == '__main__':
    unittest.main()
    print(solvers)
