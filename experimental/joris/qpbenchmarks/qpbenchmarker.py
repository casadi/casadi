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
from numpy import *
from casadi import *
from casadi.tools import *

class OnlineQPBenchMark:
  def __init__(self,name):
    self.name = name
    
    self.nQP,self.nV,self.nC,self.nEC = self.readmatrix('dims.oqp')

    self.H  = DMatrix(self.readmatrix('H.oqp'))
    self.g  = DMatrix(self.readmatrix('g.oqp'))
    self.lb = DMatrix(self.readmatrix('lb.oqp'))
    self.ub = DMatrix(self.readmatrix('ub.oqp'))

    if self.nC > 0:
        self.A   = DMatrix(self.readmatrix('A.oqp'))
        self.lbA = DMatrix(self.readmatrix('lbA.oqp'))
        self.ubA = DMatrix(self.readmatrix('ubA.oqp'))

    self.x_opt   = DMatrix(self.readmatrix('x_opt.oqp'))
    self.y_opt   = self.readmatrix('y_opt.oqp')
    self.obj_opt = self.readmatrix('obj_opt.oqp')

  def readmatrix(self,name):
    return loadtxt(self.name + '/'+name)

qp = OnlineQPBenchMark('diesel')

qpsolvers = []
try:
  qpsolvers.append(IpoptQpSolver)
except:
  pass
try:
  qpsolvers.append(OOQpSolver)
except:
  pass
try:
  qpsolvers.append(QPOasesSolver)
except:
  pass

for qpsolver in qpsolvers:
  print qpsolver

  solver = qpsolver(qp.H.sparsity(),qp.A.sparsity())
  solver.init()
  solver.setInput(qp.H,"h")
  solver.setInput(qp.A,"a")
  for i in range(qp.g.shape[1]):
    solver.setInput(qp.lbA[i,:].T,"lba")
    solver.setInput(qp.ubA[i,:].T,"uba")
    solver.setInput(qp.lb[i,:].T,"lbx")
    solver.setInput(qp.ub[i,:].T,"ubx")
    solver.setInput(qp.g[i,:].T,"g")

    solver.solve()

    print solver.getOutput("primal")
    print qp.x_opt[i,:].T
    print qp.y_opt[i,:].T
    print qp.obj_opt[i]
    print solver.getOutput("cost")
    print fabs(solver.getOutput("primal")-qp.x_opt[i,:].T)
    assert(all(fabs(solver.getOutput("primal")-qp.x_opt[i,:].T)<1e-4))
    assert(fabs(qp.obj_opt[i]-solver.getOutput("cost"))<1e-5)
