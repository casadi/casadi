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
import numpy as np
from numpy import inf, pi
import casadi as c
import numpy
import unittest
from types import *
from helpers import *

solvers= []
try:
  ca.load_rootfinder("kinsol")
  solvers.append(("kinsol",{"abstol":1e-10},[]))
except:
  pass
try:
  ca.load_nlpsol("ipopt")
  solvers.append(("nlpsol",{"nlpsol": "ipopt","nlpsol_options":{"print_time": False,"ipopt": {"print_level": 0}}},[]))
except:
  pass
try:
  solvers.append(("newton",{},[]))
except:
  pass

solvers.append(("fast_newton",{},("codegen")))

print(solvers)

class ImplicitFunctiontests(casadiTestCase):

  @memory_heavy()
  def test_linear(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      x=ca.SX.sym("x",2)
      A_ = ca.DM([[1,2],[3,2.1]])
      C_ = ca.DM([[1.6,2.1],[1,1.3]])
      b_ = ca.DM([0.7,0.6])
      f=ca.Function("f", [x],[A_ @ x-b_, C_ @ x])
      solver=ca.rootfinder("solver", Solver, f, options)
      solver_out = solver(0)

      refsol = ca.Function("refsol", [x],[ca.solve(A_,b_), C_ @ ca.solve(A_,b_)])
      self.checkfunction(solver,refsol,inputs=[0],digits=10)
      if "newton" in Solver: self.check_serialize(solver,inputs=[0])
      if "codegen" in features: self.check_codegen(solver,inputs=[0])

      A = ca.SX.sym("A",2,2)
      b = ca.SX.sym("b",2)
      f=ca.Function("f", [x,A,b],[A @ x-b])
      solver=ca.rootfinder("solver", Solver, f, options)
      solver_in = [0]*3  # type: list
      solver_in[1]=A_
      solver_in[2]=b_

      refsol = ca.Function("refsol", [x,A,b],[ca.solve(A,b)])
      self.checkfunction(solver,refsol,inputs=solver_in,digits=10)
      if "newton" in Solver: self.check_serialize(solver,inputs=solver_in)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

      A = ca.SX.sym("A",2,2)
      b = ca.SX.sym("b",2)
      f=ca.Function("f", [x,A,b],[A @ x-b,C_ @ x])
      for ad_weight_sp in [0,1]:
        for ad_weight in [0,1]:
          print(ad_weight, ad_weight_sp)
          options2 = dict(options)
          options2["ad_weight_sp"] = ad_weight_sp
          options2["ad_weight"] = ad_weight
          solver=ca.rootfinder("solver", Solver, f, options2)
          solver_in = [0]*3  # type: list
          solver_in[1]=A_
          solver_in[2]=b_

          refsol = ca.Function("refsol", [x,A,b],[ca.solve(A,b),C_ @ ca.solve(A,b)])
          self.checkfunction(solver,refsol,inputs=solver_in,digits=10)
          if "newton" in Solver: self.check_serialize(solver,inputs=solver_in)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

  def test_missing_symbols(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      x=ca.SX.sym("x")
      p=ca.SX.sym("p")
      f=ca.Function("f", [x], [ca.sin(x+p)],{"allow_free":True})
      with self.assertInException("[p] are free"):
        ca.rootfinder("solver", Solver, f, options)

  def test_scalar1(self):
    self.message("Scalar implicit problem, n=0")
    for Solver, options, features in solvers:
      self.message(Solver)
      x=ca.SX.sym("x")
      f=ca.Function("f", [x], [ca.sin(x)])
      solver=ca.rootfinder("solver", Solver, f, options)
      solver_out = solver(0)

      refsol = ca.Function("refsol", [x], [ca.ceil(x/pi-0.5)*pi])
      self.checkfunction(solver,refsol,inputs=[6],digits=10)

  def test_unsolved_stats(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      x=ca.SX.sym("x")
      f=ca.Function("f", [x], [ca.sin(x)])
      solver=ca.rootfinder("solver", Solver, f, options)
      try:
          solver.stats()
      except:
          pass # no segfault
     
  def test_scalar2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=ca.SX.sym("x")
      y=ca.SX.sym("y")
      n=0.2
      f=ca.Function("f", [y,x], [x-ca.arcsin(y)])
      solver=ca.rootfinder("solver", Solver, f, options)
      refsol = ca.Function("refsol", [y,x], [ca.sin(x)])
      self.checkfunction(solver,refsol,inputs=[0,n],digits=6,sens_der=False,failmessage=message)

  def test_scalar2_indirect(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=ca.SX.sym("x")
      y=ca.SX.sym("y")
      n=0.2
      f=ca.Function("f", [y,x], [x-ca.arcsin(y)])
      solver=ca.rootfinder("solver", Solver, f, options)

      X = ca.MX.sym("X")
      R = solver(ca.MX(),X)

      trial = ca.Function("trial", [X], [R])
      refsol = ca.Function("refsol", [x],[ca.sin(x)])
      self.checkfunction(trial,refsol,inputs=[n],digits=6,sens_der=False,failmessage=message)

  def test_large(self):
    for Solver, options, features in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue

      message = Solver
      N = 5
      s = ca.Sparsity.lower(N)
      x=ca.SX.sym("x",s)

      y=ca.SX.sym("y",s)
      y0 = ca.DM(ca.Sparsity.diag(N),0.1)

      f=ca.Function("f", [y.nz[:],x.nz[:]],[((((x+y0)) @ (x+y0).T-((y+y0)) @ (y+y0).T)[s]).nz[:]])
      options2 = dict(options)
      options2["constraints"] = [1]*s.nnz()
      solver=ca.rootfinder("options2", Solver, f, options2)

      X = ca.MX.sym("X",x.sparsity())
      R = solver(ca.MX(),X.nz[:])

      trial = ca.Function("trial", [X],[R])
      trial_in = ca.DM(trial.sparsity_in(0),[abs(ca.cos(i)) for i in range(x.nnz())])
      trial_out = trial(trial_in)

      f_in = [trial_out, trial_in.nz[:]]
      f_out = f(*f_in)

      f_in = [trial_in.nz[:], trial_in.nz[:]]
      f_out = f(*f_in)

      refsol = ca.Function("refsol", [X],[X.nz[:]])
      refsol_in = [0]*refsol.n_in()  # type: list

      refsol_in[0]=trial_in[0]

      self.checkfunction(trial,refsol,inputs=refsol_in,digits=6,sens_der=False,evals=1,failmessage=message)

  @known_bug()
  def test_vector2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=ca.SX.sym("x")
      y=ca.SX.sym("y",2)
      y0 = ca.DM([0.1,0.4])
      yy = y + y0
      n=0.2
      f=ca.Function("f", [y,x],[ca.vertcat(x-ca.arcsin(yy[0]),yy[1]**2-yy[0])])
      solver=ca.rootfinder("solver", Solver, f, options)

      refsol = ca.Function("refsol", [y,x],[ca.vertcat(ca.sin(x),ca.sqrt(ca.sin(x)))-y0]) # ,sin(x)**2])
      self.checkfunction(solver,refsol,inputs=[n,0],digits=4,sens_der=False,failmessage=message)

  def testKINSol1c(self):
    self.message("Scalar KINSol problem, n=0, constraint")
    x=ca.SX.sym("x")
    f=ca.Function("f", [x],[ca.sin(x)])
    solver=ca.rootfinder("solver", "kinsol", f, {"constraints":[-1]})
    solver_out = solver(-6)
    self.assertAlmostEqual(solver_out[0],-2*pi,5)

  def test_constraints(self):
    for Solver, options, features in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue

      print(Solver, options)
      x=ca.SX.sym("x",2)
      f=ca.Function("f", [x],[ca.vertcat(*[(x+3).T @ ((x-2)),(x-4).T @ ((x+ca.vertcat(*[1,2])))])])
      options2 = dict(options)
      options2["constraints"] = [-1,0]
      solver=ca.rootfinder("solver", Solver, f, options2)
      solver_out = solver(0)

      self.checkarray(solver_out,ca.DM([-3.0/50*(ca.sqrt(1201)-1),2.0/25*(ca.sqrt(1201)-1)]),digits=6)

      f=ca.Function("f", [x],[ca.vertcat(*[(x+3).T @ ((x-2)),(x-4).T @ ((x+ca.vertcat(*[1,2])))])])
      options2 = dict(options)
      options2["constraints"] = [1,0]
      solver=ca.rootfinder("solver", Solver, f, options2)
      solver_out = solver(0)

      self.checkarray(solver_out,ca.DM([3.0/50*(ca.sqrt(1201)+1),-2.0/25*(ca.sqrt(1201)+1)]),digits=6)

  def test_implicitbug(self):
    # Total number of variables for one finite element
    X0 = ca.MX.sym("X0")
    V = ca.MX.sym("V")

    V_eq = ca.vertcat(*[V[0]-X0])

    # Root-finding function, implicitly defines V as a function of X0 and P
    vfcn = ca.Function("vfcn", [V,X0], [V_eq], {"ad_weight":0, "ad_weight_sp":1})

    # Convert to SX to decrease overhead
    vfcn_sx = vfcn.expand('vfcn_sx', {"ad_weight":0, "ad_weight_sp":1})

    # Create a implicit function instance to solve the system of equations
    ifcn = ca.rootfinder("ifcn", "newton", vfcn_sx, {"linear_solver":"csparse"})

    #ifcn = Function('I', [X0],[vertcat(*[X0])])
    [V] = ifcn.call([0,X0],True)

    f = 1  # fails

    F = ca.Function("F", [X0], [f*X0+V], {"ad_weight":0, "ad_weight_sp":1})

    # Test values
    x0_val  = 1

    J = jacobian_old(F, 0, 0)
    J_out = J(x0_val)
    print(J_out[0])
    print(J)

    self.checkarray(J_out[0],ca.DM([2]))

  def test_extra_outputs(self):
    x = ca.SX.sym("x")
    a = ca.SX.sym("a")
    f = ca.Function("f", [x,a],[ca.tan(x)-a,ca.sqrt(a)*x**2 ])
    for Solver, options, features in solvers:
      print(Solver)
      options2 = dict(options)
      options2["ad_weight_sp"] = 1
      solver=ca.rootfinder("solver", Solver, f, options2)
      solver_in = [0.1,0.3]

      refsol = ca.Function("refsol", [x,a],[ca.arctan(a),ca.sqrt(a)*ca.arctan(a)**2])
      self.checkfunction(solver,refsol,inputs=solver_in,digits=5)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

    x = ca.SX.sym("x",2)
    a = ca.SX.sym("a",2)
    f = ca.Function("f", [x,a],[ca.tan(x)-a,ca.sqrt(a)*x**2 ])

  def test_no_success(self):

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    for Solver, options, features in solvers:
      opts = dict(options)
      opts["error_on_fail"] = False
      solver = ca.rootfinder("solver",Solver,{'x':ca.vertcat(x,y), 'g':ca.vertcat(ca.sin(x)-2,ca.sin(y)-2)},opts)
      solver(x0=0)
      self.assertFalse(solver.stats()["success"])

      opts = dict(options)
      solver = ca.rootfinder("solver",Solver,{'x':ca.vertcat(x,y), 'g':ca.vertcat(ca.sin(x)-2,ca.sin(y)-2)},opts)
      with self.assertInException("process"):
        solver(x0=0)

  def test_loop(self):
    x=ca.SX.sym("x")
    for Solver, options, features in solvers:
      solver = ca.rootfinder("solver",Solver,{"x":x,"g":x**3-2*x+2},options)
      if Solver=="kinsol": continue
      if Solver=="nlpsol": continue
      if Solver=="fast_newton": continue
      print(Solver)
      print(solver)
      res = solver(x0=0)["x"]
      self.checkarray(res,-1.7692923542386)

  def test_segfault_codegen(self):
    # Symbols
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    # B-Spline interpolant
    u = np.linspace(1, 5, 10)
    v = u ** 2
    f = ca.interpolant("interp", "bspline", [u], v, {})

    # Extrapolated interpolant, suitable for rootfinding later on
    f = ca.Function(
        "f",
        [x],
        [ca.if_else(x < u[0], x ** 2, ca.if_else(x > u[-1], x ** 2, f(x), False), False)],
    )

    # Use rootfinder to find x such that x**2=y
    res = ca.Function("res", [x, y], [y - f(x)])
    rf = ca.rootfinder("rf", "fast_newton", res)

    g = ca.Function("g",[x],[x - rf(x, x)])
    self.check_codegen(g,[1.01])
    
    
  def test_generic(self):
    a = ca.SX.sym("a")
    b = ca.SX.sym("b")
    c = ca.SX.sym("c")
    
    f = ca.Function("f",[a,b,c],[a**2-b,b**2-c,a*b*c+1],["a","b","c"],["f","g","h"])

    a0 = 1
    b0 = 2
    c0 = 3
     
    rf = ca.rootfinder("rf","newton",f,{"implicit_input":0,"implicit_output":0})
     
    res = rf(a0=a0,b=b0,c=c0)
    self.assertAlmostEqual(res["a"]**2-b0,0,10)
    self.assertAlmostEqual(b0**2-c0-res["g"],0,10)
    self.assertAlmostEqual(res["a"]*b0*c0+1-res["h"],0,10)
    print(rf)
    jac_g_x = rf.get_function("jac_g_x")
    self.assertTrue(jac_g_x.name_in()==["a","b","c"])
    self.assertTrue(jac_g_x.name_out()==["jac_f_a","f","g","h"])
    rf = ca.rootfinder("rf","newton",f,{"implicit_input":1,"implicit_output":0})
    res = rf(a=a0,b0=b0,c=c0)
    self.assertAlmostEqual(a0**2-res["b"],0,10)
    self.assertAlmostEqual(res["b"]**2-c0-res["g"],0,10)
    self.assertAlmostEqual(a0*res["b"]*c0+1-res["h"],0,10)

    jac_g_x = rf.get_function("jac_g_x")
    self.assertTrue(jac_g_x.name_in()==["a","b","c"])
    self.assertTrue(jac_g_x.name_out()==["jac_f_b","f","g","h"])
    
    rf = ca.rootfinder("rf","newton",f,{"implicit_input":2,"implicit_output":1})
    res = rf(a=a0,b=b0,c0=c0)
    self.assertAlmostEqual(a0**2-b0-res["f"],0,10)
    self.assertAlmostEqual(b0**2-res["c"],0,10)
    self.assertAlmostEqual(a0*b0*res["c"]+1-res["h"],0,10)
    jac_g_x = rf.get_function("jac_g_x")
    self.assertTrue(jac_g_x.name_in()==["a","b","c"])
    self.assertTrue(jac_g_x.name_out()==["jac_g_c","f","g","h"])

    rf = ca.rootfinder("rf","newton",{"x":a,"p":ca.vertcat(b,c),"g":a**2-b})
    res = rf(x0=a0,p=ca.vertcat(b0,c0))
    self.assertAlmostEqual(res["x"]**2-b0,0,10)
    jac_g_x = rf.get_function("jac_g_x")
    self.assertTrue(jac_g_x.name_in()==["x","p"])
    self.assertTrue(jac_g_x.name_out()==["jac_g_x","g"])
    
if __name__ == '__main__':
    unittest.main()
