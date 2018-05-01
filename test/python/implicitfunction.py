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

solvers= []
try:
  load_rootfinder("kinsol")
  solvers.append(("kinsol",{"abstol":1e-10},[]))
except:
  pass
try:
  load_nlpsol("ipopt")
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
      x=SX.sym("x",2)
      A_ = DM([[1,2],[3,2.1]])
      C_ = DM([[1.6,2.1],[1,1.3]])
      b_ = DM([0.7,0.6])
      f=Function("f", [x],[mtimes(A_,x)-b_, mtimes(C_,x)])
      solver=rootfinder("solver", Solver, f, options)
      solver_out = solver(0)

      refsol = Function("refsol", [x],[solve(A_,b_), mtimes(C_,solve(A_,b_))])
      self.checkfunction(solver,refsol,inputs=[0],digits=10)
      if "codegen" in features: self.check_codegen(solver,inputs=[0])

      A = SX.sym("A",2,2)
      b = SX.sym("b",2)
      f=Function("f", [x,A,b],[mtimes(A,x)-b])
      solver=rootfinder("solver", Solver, f, options)
      solver_in = [0]*3
      solver_in[1]=A_
      solver_in[2]=b_

      refsol = Function("refsol", [x,A,b],[solve(A,b)])
      self.checkfunction(solver,refsol,inputs=solver_in,digits=10)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

      A = SX.sym("A",2,2)
      b = SX.sym("b",2)
      f=Function("f", [x,A,b],[mtimes(A,x)-b,mtimes(C_,x)])
      for ad_weight_sp in [0,1]:
        for ad_weight in [0,1]:
          print(ad_weight, ad_weight_sp)
          options2 = dict(options)
          options2["ad_weight_sp"] = ad_weight_sp
          options2["ad_weight"] = ad_weight
          solver=rootfinder("solver", Solver, f, options2)
          solver_in = [0]*3
          solver_in[1]=A_
          solver_in[2]=b_

          refsol = Function("refsol", [x,A,b],[solve(A,b),mtimes(C_,solve(A,b))])
          self.checkfunction(solver,refsol,inputs=solver_in,digits=10)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

  def test_missing_symbols(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      x=SX.sym("x")
      p=SX.sym("p")
      f=Function("f", [x], [sin(x+p)])
      with self.assertInException("[p] are free"):
        rootfinder("solver", Solver, f, options)

  def test_scalar1(self):
    self.message("Scalar implicit problem, n=0")
    for Solver, options, features in solvers:
      self.message(Solver)
      x=SX.sym("x")
      f=Function("f", [x], [sin(x)])
      solver=rootfinder("solver", Solver, f, options)
      solver_out = solver(0)

      refsol = Function("refsol", [x], [ceil(x/pi-0.5)*pi])
      self.checkfunction(solver,refsol,inputs=[6],digits=10)

  def test_scalar2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=SX.sym("x")
      y=SX.sym("y")
      n=0.2
      f=Function("f", [y,x], [x-arcsin(y)])
      solver=rootfinder("solver", Solver, f, options)
      refsol = Function("refsol", [y,x], [sin(x)])
      self.checkfunction(solver,refsol,inputs=[0,n],digits=6,sens_der=False,failmessage=message)

  def test_scalar2_indirect(self):
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=SX.sym("x")
      y=SX.sym("y")
      n=0.2
      f=Function("f", [y,x], [x-arcsin(y)])
      solver=rootfinder("solver", Solver, f, options)

      X = MX.sym("X")
      R = solver(MX(),X)

      trial = Function("trial", [X], [R])
      refsol = Function("refsol", [x],[sin(x)])
      self.checkfunction(trial,refsol,inputs=[n],digits=6,sens_der=False,failmessage=message)

  def test_large(self):
    for Solver, options, features in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue

      message = Solver
      N = 5
      s = Sparsity.lower(N)
      x=SX.sym("x",s)

      y=SX.sym("y",s)
      y0 = DM(Sparsity.diag(N),0.1)

      f=Function("f", [y.nz[:],x.nz[:]],[((mtimes((x+y0),(x+y0).T)-mtimes((y+y0),(y+y0).T))[s]).nz[:]])
      options2 = dict(options)
      options2["constraints"] = [1]*s.nnz()
      solver=rootfinder("options2", Solver, f, options2)

      X = MX.sym("X",x.sparsity())
      R = solver(MX(),X.nz[:])

      trial = Function("trial", [X],[R])
      trial_in = DM(trial.sparsity_in(0),[abs(cos(i)) for i in range(x.nnz())])
      trial_out = trial(trial_in)

      f_in = [trial_out, trial_in.nz[:]]
      f_out = f(*f_in)

      f_in = [trial_in.nz[:], trial_in.nz[:]]
      f_out = f(*f_in)

      refsol = Function("refsol", [X],[X.nz[:]])
      refsol_in = [0]*refsol.n_in();refsol_in[0]=trial_in[0]

      self.checkfunction(trial,refsol,inputs=refsol_in,digits=6,sens_der=False,evals=1,failmessage=message)

  @known_bug()
  def test_vector2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options, features in solvers:
      self.message(Solver)
      message = Solver
      x=SX.sym("x")
      y=SX.sym("y",2)
      y0 = DM([0.1,0.4])
      yy = y + y0
      n=0.2
      f=Function("f", [y,x],[vertcat(*[x-arcsin(yy[0]),yy[1]**2-yy[0]])])
      solver=rootfinder("solver", Solver, f, options)
      solver_in = [0]*solver.n_in();solver_in[0]=n
      solver_out = solver(solver_in)

      refsol = Function("refsol", [y,x],[vertcat(*[sin(x),sqrt(sin(x))])-y0]) # ,sin(x)**2])
      refsol_in = [0]*refsol.n_in();refsol_in[0]=n
      self.checkfunction(solver,refsol,digits=5,sens_der=False,failmessage=message)

  def testKINSol1c(self):
    self.message("Scalar KINSol problem, n=0, constraint")
    x=SX.sym("x")
    f=Function("f", [x],[sin(x)])
    solver=rootfinder("solver", "kinsol", f, {"constraints":[-1]})
    solver_out = solver(-6)
    self.assertAlmostEqual(solver_out[0],-2*pi,5)

  def test_constraints(self):
    for Solver, options, features in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue

      print(Solver, options)
      x=SX.sym("x",2)
      f=Function("f", [x],[vertcat(*[mtimes((x+3).T,(x-2)),mtimes((x-4).T,(x+vertcat(*[1,2])))])])
      options2 = dict(options)
      options2["constraints"] = [-1,0]
      solver=rootfinder("solver", Solver, f, options2)
      solver_out = solver(0)

      self.checkarray(solver_out,DM([-3.0/50*(sqrt(1201)-1),2.0/25*(sqrt(1201)-1)]),digits=6)

      f=Function("f", [x],[vertcat(*[mtimes((x+3).T,(x-2)),mtimes((x-4).T,(x+vertcat(*[1,2])))])])
      options2 = dict(options)
      options2["constraints"] = [1,0]
      solver=rootfinder("solver", Solver, f, options2)
      solver_out = solver(0)

      self.checkarray(solver_out,DM([3.0/50*(sqrt(1201)+1),-2.0/25*(sqrt(1201)+1)]),digits=6)

  def test_implicitbug(self):
    # Total number of variables for one finite element
    X0 = MX.sym("X0")
    V = MX.sym("V")

    V_eq = vertcat(*[V[0]-X0])

    # Root-finding function, implicitly defines V as a function of X0 and P
    vfcn = Function("vfcn", [V,X0], [V_eq], {"ad_weight":0, "ad_weight_sp":1})

    # Convert to SX to decrease overhead
    vfcn_sx = vfcn.expand('vfcn_sx', {"ad_weight":0, "ad_weight_sp":1})

    # Create a implicit function instance to solve the system of equations
    ifcn = rootfinder("ifcn", "newton", vfcn_sx, {"linear_solver":"csparse"})

    #ifcn = Function('I', [X0],[vertcat(*[X0])])
    [V] = ifcn.call([0,X0],True)

    f = 1  # fails

    F = Function("F", [X0], [f*X0+V], {"ad_weight":0, "ad_weight_sp":1})

    # Test values
    x0_val  = 1

    J = F.jacobian_old(0, 0)
    J_out = J(x0_val)
    print(J_out[0])
    print(J)

    self.checkarray(J_out[0],DM([2]))

  def test_extra_outputs(self):
    x = SX.sym("x")
    a = SX.sym("a")
    f = Function("f", [x,a],[tan(x)-a,sqrt(a)*x**2 ])
    for Solver, options, features in solvers:
      print(Solver)
      options2 = dict(options)
      options2["ad_weight_sp"] = 1
      solver=rootfinder("solver", Solver, f, options2)
      solver_in = [0.1,0.3]

      refsol = Function("refsol", [x,a],[arctan(a),sqrt(a)*arctan(a)**2])
      self.checkfunction(solver,refsol,inputs=solver_in,digits=5)
      if "codegen" in features: self.check_codegen(solver,inputs=solver_in)

    x = SX.sym("x",2)
    a = SX.sym("a",2)
    f = Function("f", [x,a],[tan(x)-a,sqrt(a)*x**2 ])

  def test_no_success(self):

    x=SX.sym("x")
    y=SX.sym("y")

    for Solver, options, features in solvers:
      opts = dict(options)
      if Solver=="kinsol": opts["error_on_fail"] = False # has different default
      solver = rootfinder("solver",Solver,{'x':vertcat(x,y), 'g':vertcat(sin(x)-2,sin(y)-2)},opts)
      solver(x0=0)
      self.assertFalse(solver.stats()["success"])

      opts = dict(options)
      opts["error_on_fail"] = True
      solver = rootfinder("solver",Solver,{'x':vertcat(x,y), 'g':vertcat(sin(x)-2,sin(y)-2)},opts)
      with self.assertInException("process"):
        solver(x0=0)

if __name__ == '__main__':
    unittest.main()
