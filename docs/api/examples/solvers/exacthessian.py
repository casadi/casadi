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
#! Exact Hessian
#! =====================
from casadi import *
from numpy import *
import casadi as c

#! We will investigate the use of an exact Hessian with the help of the Rosenbrock function
x=SX.sym("x")
y=SX.sym("y")
    
obj = (1-x)**2+100*(y-x**2)**2
#! We choose to add a single constraint
constr = x**2+y**2

nlp=SXFunction(nlpIn(x=vertcat([x,y])),nlpOut(f=obj,g=constr))
solver = NlpSolver("ipopt", nlp)
    
#! We need the hessian of the lagrangian.
#! A problem with n decision variables and m constraints gives us a hessian of size n x n
  
sigma=SX.sym("sigma")  # A scalar factor
lambd=SX.sym("lambd")  # Multipier of the problem, shape m x 1.

xy = vertcat([x,y])

h=SXFunction(hessLagIn(x=xy,lam_g=lambd,lam_f=sigma),
             hessLagOut(hess=sigma*hessian(obj,xy)+lambd*hessian(constr,xy)))
   
#! We solve the problem with an exact hessian
solver = NlpSolver("ipopt", nlp)
solver.setOption("hess_lag",h)
solver.init()
solver.setInput([-10]*2,"lbx")
solver.setInput([10]*2,"ubx")
solver.setInput([0],"lbg")
solver.setInput([1],"ubg")
solver.evaluate()

for sol in array(solver.getOutput()):
  print "%.15f" % sol

#! To compare the behaviour of convergence, we solve the same problem without exact hessian
solver = NlpSolver("ipopt", nlp)
solver.init()
solver.setInput([-10]*2,"lbx")
solver.setInput([10]*2,"ubx")
solver.setInput([0],"lbg")
solver.setInput([1],"ubg")
solver.evaluate()

for sol in array(solver.getOutput()):
  print "%.15f" % sol

