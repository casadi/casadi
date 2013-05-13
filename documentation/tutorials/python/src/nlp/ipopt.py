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
#! CasADi tutorial
#! ==================
#! This tutorial file explains the interface with IPOPT
#! Ipopt solves problems of the form:
#!
#!
#! Minimize     f(x)
#! x in R^n
#! s.t         g_L <= g(x) <= g_U
#! x_L <= x    <= x_U
from numpy import *
import numpy as n
from casadi import *
#! Let's solve a simple scalar non-linear program:

x = ssym("x")

y = ssym("y")

#f  = SXFunction([x,y], tan(x)-1) 

#g=SXFunction([x,y],[y])

#print f.eval(x)
#! Quadractic program
#! ------------------
#! (Ab)using Ipopt to do a simple quadratic problem
#!
#! 
P = n.eye(5)
A = n.diag([2.0,1,4,10,-2])
q = [1,0,0,0,0]
b = [1,1,1,1,1]

X = MX("x",5,1)
P = MX(DMatrix(P))
q = MX(DMatrix(q))
A = MX(DMatrix(A))

#! Objective
F = 0.5*mul(mul(trans(X),P),X) + mul(trans(q),X)

#! Constraint
G = X+X

#! NLP
nlp = MXFunction(nlpIn(x=X),nlpOut(f=F,g=G))
nlp.init()
nlp.setInput([1,1,1,1,1],"x")
nlp.evaluate()
#! Test the objective for some value of x:
print nlp.output("f").toArray()

solver = IpoptSolver(nlp)
solver.printOptions()

solver.init()

#! The default lower an upper bound on the optimizations variables is zero.
#! Change them to unbouded as follows:
solver.setInput([-100,-100,-100,-100,-100],"lbx")
solver.setInput([100, 100, 100, 100, 100],"ubx")

#! Inequality constraints.
#! The lower bound is also necessary, although Ipopt itself does not seem to require it
solver.setInput(b,"ubg")
solver.setInput([-100,-100,-100,-100,-100],"lbg")


solver.solve()
print solver.output("x")
#! Nested optimization
