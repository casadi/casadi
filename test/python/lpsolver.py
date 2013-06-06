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

lpsolvers = []
try:
  lpsolvers.append((QPLPSolver,{"qp_solver":OOQPSolver,"qp_solver_options": {"mutol": 1e-12, "artol": 1e-12} }))
except:
  pass
  
print lpsolvers

class LPSolverTests(casadiTestCase):

  def test_all(self):

    A = DMatrix([ [-1,1],[1,1],[1,-2]])
    LBA = DMatrix([ -inf, 2, -inf ])
    UBA = DMatrix([ 1, inf, 4 ])
    LBX = DMatrix([ -inf, 0 ])
    UBX = DMatrix([ inf, inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for lpsolver, lp_options in lpsolvers:
      self.message("lpsolver: " + str(lpsolver))

      solver = lpsolver(A.sparsity())
      solver.setOption(lp_options)
      solver.init()

      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.solve()

      self.checkarray(solver.getOutput(),DMatrix([0.5,1.5]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0.5,-1.5,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],2.5,5,str(lpsolver))
      
      # Make a LBX active
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput([-inf,2],"lbx")
      solver.setInput(UBX,"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.solve()

      self.checkarray(solver.getOutput(),DMatrix([1,2]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,-3]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([2,0,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],4,5,str(lpsolver))
      
      # Make a UBX active
      solver.setInput(c,"c")
      solver.setInput(A,"a")
      solver.setInput(LBX,"lbx")
      solver.setInput([inf,1],"ubx")
      solver.setInput(LBA,"lba")
      solver.setInput(UBA,"uba")

      solver.solve()

      self.checkarray(solver.getOutput(),DMatrix([1,1]),str(lpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,1]),str(lpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([0,-2,0]),str(lpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],3,5,str(lpsolver))
      
if __name__ == '__main__':
    unittest.main()
