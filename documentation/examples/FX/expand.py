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
#! expand
#!======================
from casadi import *
import casadi as c

#! We construct a simple MX expression
x = MX.sym("x",2,2)
y = MX.sym("y",2,1)

z = mul(x,y)

#! Let's construct an MXfunction
f = MXFunction([x,y],[z])
f.init()

#! An MX graph is lazy in evaluation
print "Expression = ", f.outputExpr(0)

#! We expand the MXFunction into an SXFunction
fSX = f.expand()

print "Expanded expression = ", fSX.outputExpr(0)


#! Limitations
#! =============
#! Not all MX graphs can be expanded.
#! Here is an example of a situation where it will not work.
#!
linear_solver = CSparse(x.sparsity())
linear_solver.init()
g = linear_solver.solve(x,y)
G = MXFunction([x,y],[g])
G.init()

#! This function cannot be expanded.
try:
  G.expand()
except Exception as e:
  print e

