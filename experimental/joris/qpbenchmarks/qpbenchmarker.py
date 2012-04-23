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
from numpy import *
from casadi import *
from casadi.tools import *

class OnlineQPBenchMark:
  def __init__(self,name):
    self.name = name
    
    self.nQP,self.nV,self.nC,self.nEC = self.readmatrix('dims.oqp')

    self.H  = self.readmatrix('H.oqp')
    self.g  = self.readmatrix('g.oqp')
    self.lb = self.readmatrix('lb.oqp')
    self.ub = self.readmatrix('ub.oqp')

    if self.nC > 0:
        self.A   = self.readmatrix('A.oqp')
        self.lbA = self.readmatrix('lbA.oqp')
        self.ubA = self.readmatrix('ubA.oqp')

    self.x_opt   = self.readmatrix('x_opt.oqp')
    self.y_opt   = self.readmatrix('y_opt.oqp')
    self.obj_opt = self.readmatrix('obj_opt.oqp')

  def readmatrix(self,name):
    return loadtxt(self.name + '/'+name)

qp = OnlineQPBenchMark('diesel')

V = Variables()

V.x = ssym("x",int(qp.nV))
#V.g = ssym("g",int(qp.nV))

V.freeze()

x = V.x
g = SXMatrix(qp.g[0,:])
H = qp.H
A = qp.A

f = SXFunction([V.veccat()],[0.5*mul(mul(x.T,H),x) + mul(g.T,x)])
g = SXFunction([V.veccat()],[mul(A,x)])
g.init()

solver = WorhpSolver(f,g)
#solver.setOption("tol",1e-12)
#solver.setOption("hessian_approximation","exact")
solver.setOption("generate_hessian",True)
solver.setOption("UserHM",True)
solver.setOption("TolOpti",1e-12)
#solver.setOption("monitor",["eval_h","eval_f","eval_g","eval_jac_g","eval_grad_f"])
solver.init()

for i in range(1):
  solver.input(NLP_LBX)[V.i_x] = qp.lb[i,:].T
  solver.input(NLP_UBX)[V.i_x] = qp.ub[i,:].T
  solver.input(NLP_LBG).set(qp.lbA[i,:].T)
  solver.input(NLP_UBG).set(qp.ubA[i,:].T)
  #solver.input(NLP_LBX)[V.i_g] = qp.g[i,:].T
  #solver.input(NLP_UBX)[V.i_g] = qp.g[i,:].T
  #solver.input(NLP_X_INIT)[V.i_g] = qp.g[i,:].T
  solver.evaluate()
  print "Error = ", max(fabs(solver.output()[V.i_x] - qp.x_opt[i,:]))
  
