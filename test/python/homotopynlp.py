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
import itertools

#CasadiOptions.setCatchErrorsPython(False)

solvers= []
 
try:
  solvers.append((SimpleHomotopyNlpSolver,{"nlp_solver":"ipopt","nlp_solver_options": {"tol": 1e-12} }))
  print "Will test SimpleHomotopyNlpSolver"
except:
  pass
  

class NLPtests(casadiTestCase):

  def test_simple(self):
    x = SX.sym("x")
    y = SX.sym("y")
    tau = SX.sym("tau")

    nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=(1-x)**2+100*(y-x**2)**2))

    hnlp = SXFunction(hnlpIn(x=vertcat([x,y]),tau=tau),nlpOut(f=(1-x)**2+100*(y-x**2)**2+x-tau,g=x**2+y**2 - (tau -1.2)**2))
    hnlp.init()

    for Solver, solver_options in solvers:
      hnlpsolver = Solver(hnlp)

      hnlpsolver.setOption(solver_options)

      hnlpsolver.init()

      hnlpsolver.setInput(-Inf,"lbg")
      hnlpsolver.setInput(0,"ubg")

      hnlpsolver.evaluate()

      self.checkarray(hnlpsolver.getOutput("x"),DMatrix([1.9635508794099052e-01,3.8009070491114780e-02]))
      self.checkarray(hnlpsolver.getOutput("lam_x"),DMatrix([0,0]))
      self.checkarray(hnlpsolver.getOutput("lam_g"),DMatrix([1.4371571368129470e+00]))
      
if __name__ == '__main__':
    unittest.main()
    print solvers

