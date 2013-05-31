#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
  solver.append((KinsolSolver,{"linear_solver": CSparse}))
except:
  pass
try:
  solver.append((NLPImplicitSolver,{"linear_solver": CSparse,"nlp_solver": IpoptSolver}))
except:
  pass
try:
  solver.append((NewtonImplicitSolver,{"linear_solver": CSparse}))
except:
  pass

class NLPtests(casadiTestCase):
  def test_scalar1(self):
    self.message("Scalar implicit problem, n=0")
    for Solver, options in solvers:
      self.message(Solver.__name__)
      x=SX("x")
      f=SXFunction([x],[sin(x)])
      f.init()
      solver=Solver(f)
      solver.setOption(options)
      solver.init()
      solver.setOutput(6)
      solver.solve()
      
      refsol = SXFunction([],[2*pi])
      refsol.init()
      self.checkfx(solver,refsol,digits=5,gradient=False,hessian=False,sens_der=False)         
      
  def test_scalar2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options in solvers:
      self.message(Solver.__name__)
      message = Solver.__name__
      x=SX("x")
      y=SX("y")
      n=0.2
      f=SXFunction([y,x],[x-arcsin(y)])
      f.init()
      solver=Solver(f)
      solver.setOption(options)
      solver.init()
      solver.setFwdSeed(1)
      solver.setAdjSeed(1)
      solver.setInput(n)
      solver.evaluate(1,1)
      
      refsol = SXFunction([x],[sin(x)])
      refsol.init()
      refsol.setInput(n)
      self.checkfx(solver,refsol,digits=6,gradient=False,hessian=False,sens_der=False,failmessage=message)
      
      
  def test_vector2(self):
    self.message("Scalar implicit problem, n=1")
    for Solver, options in solvers:
      self.message(Solver.__name__)
      message = Solver.__name__
      x=SX("x")
      y=ssym("y",2)
      n=0.2
      f=SXFunction([y,x],[vertcat([x-arcsin(y[0]),y[1]**2-y[0]])])
      f.init()
      solver=Solver(f)
      solver.setOption(options)
      solver.init()
      solver.setFwdSeed(1)
      solver.setAdjSeed(1)
      solver.setInput(n)
      solver.setOutput([0.1,0.4])
      solver.evaluate(1,1)
      
      refsol = SXFunction([x],[vertcat([sin(x),sqrt(sin(x))])]) # ,sin(x)**2])
      refsol.init()
      refsol.setInput(n)
      self.checkfx(solver,refsol,digits=6,gradient=False,hessian=False,sens_der=False,failmessage=message)
      
  def testKINSol1c(self):
    self.message("Scalar KINSol problem, n=0, constraint")
    x=SX("x")
    f=SXFunction([x],[sin(x)])
    f.init()
    solver=KinsolSolver(f)
    solver.setOption("constraints",[-1])
    print solver.dictionary()
    solver.init()
    solver.setOutput(-6)
    solver.solve()
    self.assertAlmostEqual(solver.output()[0],-2*pi,5)
    
if __name__ == '__main__':
    unittest.main()

