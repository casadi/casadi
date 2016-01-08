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

#! Quadractic program
#! ------------------
#! Using Ipopt to do a simple quadratic problem
#!
#! 
P = n.eye(5)
A = n.diag([2.0,1,4,10,-2])
q = [1,0,0,0,0]
b = [1,1,1,1,1]

X = MX.sym("x",5,1)
P = MX(DM(P))
q = MX(DM(q))
A = MX(DM(A))

#! Objective
F = 0.5*mtimes([X.T,P,X]) + mtimes(q.T,X)

#! Constraint
G = X+X

#! NLP
nlp = {'x':X, 'f':F, 'g':G}

solver = nlpsol("nlp","ipopt", nlp)
solver.printOptions()

#! The default lower an upper bound on the optimizations variables is zero.
#! Change them to unbouded as follows:
lbx = [-100,-100,-100,-100,-100]
ubx = [ 100, 100, 100, 100, 100]

#! Inequality constraints.
#! The lower bound is also necessary, although Ipopt itself does not seem to require it
ubg = b
lbg = [-100,-100,-100,-100,-100]

sol = solver({"lbx":lbx, "ubx":ubx, "lbg":lbg, "ubg":ubg})
print sol["x"]
#! Nested optimization
