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
import numpy as NP

x = ssym("x")
y = ssym("y")
z = ssym("z")
v = vertcat([x,y,z])

f = SXFunction([v],[x**2 + 100*z**2])
g = SXFunction([v],[z + (1-x)**2 - y])

solv = IpoptSolver(f,g)

class Log:
  def __init__(self):
    self.iter = 0 
  def __call__(self,f,*args):
    print "====Hey, I'm an iteration===="
    print "X_OPT = ", f.getInput("x")
    print f.getStats()
    self.iter = self.iter + 1
    if self.iter > 5:
      print "This is quite enough."
      f.setOutput(1,0)

log = Log()

c = PyFunction( log, [sp_dense(3,1),sp_dense(1,1),sp_dense(1,1),sp_dense(3,1)], [sp_dense(1,1)] )
c.init()

solv.setOption("iteration_callback",c)

solv.init()

solv.setInput([2.5,3.0,0.75],"x0")
solv.setInput(0,"ubg")
solv.setInput(0,"lbg")
solv.solve()

print solv.getOutput("x")

