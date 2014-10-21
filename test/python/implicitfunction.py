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

solvers= []
try:
  LinearSolver.loadPlugin("csparse")
  ImplicitFunction.loadPlugin("kinsol")
  solvers.append(("kinsol",{"linear_solver": "csparse","abstol":1e-10}))
except:
  pass
try:
  LinearSolver.loadPlugin("csparse")
  NlpSolver.loadPlugin("ipopt")
  solvers.append(("nlp",{"linear_solver": "csparse", "nlp_solver": "ipopt"}))
except:
  pass
try:
  LinearSolver.loadPlugin("csparse")
  solvers.append(("newton",{"linear_solver": "csparse"}))
except:
  pass

class NLPtests(casadiTestCase):
  def test_scalar1(self):
    self.message("Scalar implicit problem, n=0")
    for Solver, options in solvers:
      self.message(Solver)
      x=SX.sym("x")
      f=SXFunction([x],[sin(x)])
      f.init()
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.init()
      solver.setInput(6)
      solver.evaluate()
      
      refsol = SXFunction([x],[ceil(x/pi-0.5)*pi])
      refsol.init()
      refsol.setInput(6)
      self.checkfunction(solver,refsol,digits=5)         
      
  def test_scalar2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options in solvers:
      self.message(Solver)
      message = Solver
      x=SX.sym("x")
      y=SX.sym("y")
      n=0.2
      f=SXFunction([y,x],[x-arcsin(y)])
      f.init()
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.init()
      solver.setInput(n,1)

      refsol = SXFunction([y,x],[sin(x)])
      refsol.init()
      refsol.setInput(n,1)
      self.checkfunction(solver,refsol,digits=6,sens_der=False,failmessage=message)

  def test_scalar2_indirect(self):
    for Solver, options in solvers:
      self.message(Solver)
      message = Solver
      x=SXElement.sym("x")
      y=SXElement.sym("y")
      n=0.2
      f=SXFunction([y,x],[x-arcsin(y)])
      f.init()
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.init()
      
      X = MX.sym("X")
      [R] = solver.call([MX(),X])
      
      trial = MXFunction([X],[R])
      trial.init()
      trial.setInput(n)
      
      refsol = SXFunction([x],[sin(x)])
      refsol.init()
      refsol.setInput(n)
      self.checkfunction(trial,refsol,digits=6,sens_der=False,failmessage=message)
      
  def test_large(self):
    for Solver, options in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue
      
      message = Solver
      N = 5
      s = Sparsity.tril(N)
      x=SX.sym("x",s)

      y=SX.sym("y",s)
      y0 = DMatrix(Sparsity.diag(N),0.1)

      f=SXFunction([vecNZ(y),vecNZ(x)],[vecNZ((mul((x+y0),(x+y0).T)-mul((y+y0),(y+y0).T))[s])])
      f.init()
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      # Cholesky is only unique for positive diagonal entries
      solver.setOption("constraints",[1]*s.size())
      solver.init()
      
      X = MX.sym("X",x.sparsity())
      [R] = solver.call([MX(),vecNZ(X)])
      
      trial = MXFunction([X],[R])
      trial.init()
      trial.setInput([abs(cos(i)) for i in range(x.size())])
      trial.evaluate()

      f.setInput(trial.getOutput(),0)
      f.setInput(vecNZ(trial.getInput()),1)
      f.evaluate()

      f.setInput(vecNZ(trial.getInput()),0)
      f.setInput(vecNZ(trial.getInput()),1)
      f.evaluate()
      
      refsol = MXFunction([X],[vecNZ(X)])
      refsol.init()
      refsol.setInput(trial.getInput())

      self.checkfunction(trial,refsol,digits=6,sens_der=False,evals=1,failmessage=message)
      
  @known_bug()  
  def test_vector2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options in solvers:
      self.message(Solver)
      message = Solver
      x=SXElement.sym("x")
      y=SX.sym("y",2)
      y0 = DMatrix([0.1,0.4])
      yy = y + y0
      n=0.2
      f=SXFunction([y,x],[vertcat([x-arcsin(yy[0]),yy[1]**2-yy[0]])])
      f.init()
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.init()
      solver.setInput(n)
      solver.evaluate()
      
      refsol = SXFunction([y,x],[vertcat([sin(x),sqrt(sin(x))])-y0]) # ,sin(x)**2])
      refsol.init()
      refsol.setInput(n)
      self.checkfunction(solver,refsol,digits=5,sens_der=False,failmessage=message)
      
  def testKINSol1c(self):
    self.message("Scalar KINSol problem, n=0, constraint")
    x=SXElement.sym("x")
    f=SXFunction([x],[sin(x)])
    f.init()
    solver=ImplicitFunction("kinsol",f)
    solver.setOption("constraints",[-1])
    #print solver.dictionary()
    solver.init()
    solver.setInput(-6)
    solver.evaluate()
    self.assertAlmostEqual(solver.getOutput()[0],-2*pi,5)
    
  def test_constraints(self):
    for Solver, options in solvers:
      if 'kinsol' in str(Solver): continue
      if 'newton' in str(Solver): continue
      x=SX.sym("x",2)
      f=SXFunction([x],[vertcat([mul((x+3).T,(x-2)),mul((x-4).T,(x+vertcat([1,2])))])])
      f.init()
      
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.setOption("constraints",[-1,0])
      solver.init()
      solver.evaluate()
      
      self.checkarray(solver.getOutput(),DMatrix([-3.0/50*(sqrt(1201)-1),2.0/25*(sqrt(1201)-1)]),digits=6)

      f=SXFunction([x],[vertcat([mul((x+3).T,(x-2)),mul((x-4).T,(x+vertcat([1,2])))])])
      f.init()
      
      solver=ImplicitFunction(Solver,f)
      solver.setOption(options)
      solver.setOption("constraints",[1,0])
      solver.init()
      solver.evaluate()
      
      self.checkarray(solver.getOutput(),DMatrix([3.0/50*(sqrt(1201)+1),-2.0/25*(sqrt(1201)+1)]),digits=6)

  def test_implicitbug(self):
    # Total number of variables for one finite element
    X0 = MX.sym("X0")
    V = MX.sym("V")

    V_eq = vertcat([V[0]-X0])

    # Root-finding function, implicitly defines V as a function of X0 and P
    vfcn = MXFunction([V,X0],[V_eq])
    vfcn.init()

    # Convert to SXFunction to decrease overhead
    vfcn_sx = SXFunction(vfcn)
    vfcn_sx.setOption("name","S")
    vfcn_sx.setOption("ad_mode","forward")

    # Create a implicit function instance to solve the system of equations
    ifcn = ImplicitFunction("newton",vfcn_sx)
    ifcn.setOption("linear_solver","csparse")
    ifcn.init()

    #ifcn = MXFunction([X0],[vertcat([X0])])
    #ifcn.setOption("name","I")
    #ifcn.init()
    [V] = ifcn.call([0,X0],True)

    f = 1  # fails

    F = MXFunction([X0],[f*X0+V])
    F.setOption("name","F")
    F.setOption("ad_mode","forward")
    F.init()

    # Test values
    x0_val  = 1

    G = F.gradient(0,0)
    G.init()
    G.setInput(x0_val)
    G.evaluate()
    print G.getOutput()
    print G

    J = F.jacobian(0,0)
    J.init()
    J.setInput(x0_val)
    J.evaluate()
    print J.getOutput()
    print J
    
    self.checkarray(G.getOutput(),DMatrix([2]))
    self.checkarray(J.getOutput(),DMatrix([2]))
    
if __name__ == '__main__':
    unittest.main()

