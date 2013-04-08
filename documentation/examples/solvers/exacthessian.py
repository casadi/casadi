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
#! Exact Hessian
#! =====================
from casadi import *
from numpy import *
import casadi as c

#! We will investigate the use of an exact Hessian with the help of the Rosenbrock function
x=SX("x")
y=SX("y")
    
obj = (1-x)**2+100*(y-x**2)**2
#! We choose to add a single constraint
constr = x**2+y**2

f=SXFunction([vertcat([x,y])],[obj])
g=SXFunction([vertcat([x,y])],[constr])
solver = IpoptSolver(f,g)
    
#! We need the hessian of the lagrangian.
#! A problem with n decision variables and m constraints gives us a hessian of size n x n
  
sigma=SX("sigma")  # A scalar factor
lambd=SX("lambd")  # Multipier of the problem, shape m x 1.

xy = vertcat([x,y])

h=SXFunction([xy,lambd,sigma],[sigma*hessian(obj,xy)+lambd*hessian(constr,xy)])
   
#! We solve the problem with an exact hessian
solver = IpoptSolver(f,g,h)
solver.init()
solver.input("lbx").set([-10]*2)
solver.input("ubx").set([10]*2)
solver.input("lbg").set([0])
solver.input("ubg").set([1])
solver.solve()

for sol in array(solver.output()):
  print "%.15f" % sol

#! To compare the behaviour of convergence, we solve the same problem without exact hessian
solver = IpoptSolver(f,g)
solver.init()
solver.input("lbx").set([-10]*2)
solver.input("ubx").set([10]*2)
solver.input("lbg").set([0])
solver.input("ubg").set([1])
solver.solve()

for sol in array(solver.output()):
  print "%.15f" % sol

