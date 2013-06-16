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

qcqpsolvers = []
try:
  qcqpsolvers.append((SOCPQCQPSolver,{"socp_solver": SDPSOCPSolver, "socp_solver_options": {"sdp_solver": DSDPSolver} },False))
except:
  pass


class QCQPSolverTests(casadiTestCase):

  def test_bounds(self):
    #  min  2 x + y
    #   x,y
    #
    #  s.t.  x^2 + 2y^2 + 2*x + 3*y - 7 <= 0
    H = DMatrix.zeros(2,2)
    G = DMatrix([2,1])
    A = DMatrix(0,2)
    P = DMatrix([[1,0],[0,1]])
    Q = DMatrix([2,3])
    R = DMatrix([-7])
    LBX = DMatrix([ -inf, -inf ])
    UBX = DMatrix([ inf, inf ])
    c = DMatrix([ 2.0, 1.0 ])
    
    for qcqpsolver, qcqp_options, re_init in qcqpsolvers:
      self.message("qcqpsolver: " + str(qcqpsolver))

      solver = qcqpsolver(qcqpStruct(a=A.sparsity(),p=P.sparsity(),h=H.sparsity()))
      solver.setOption(qcqp_options)
      solver.init()

      solver.setInput(H,"h")
      solver.setInput(G,"g")
      solver.setInput(A,"a")
      solver.setInput(P,"p")
      solver.setInput(Q,"q")
      solver.setInput(R,"r")
      solver.setInput(LBX,"lbx")
      solver.setInput(UBX,"ubx")

      solver.solve()

      return
      self.checkarray(solver.getOutput(),DMatrix([-(sqrt(73)+3)/3,-(sqrt(73)+9)/12]),str(qcqpsolver),digits=5)
      self.checkarray(solver.getOutput("lam_x"),DMatrix([0,0]),str(qcqpsolver),digits=5)

      self.checkarray(solver.getOutput("lam_a"),DMatrix([]),str(qcqpsolver),digits=5)
      
      self.assertAlmostEqual(solver.getOutput("cost")[0],0,5,str(qcqpsolver))
      
if __name__ == '__main__':
    unittest.main()
