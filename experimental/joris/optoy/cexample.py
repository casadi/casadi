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
from coptoy import *
import numpy as np

# Example 1: SDP problem
x = var()
y = var()

A = blockcat([[-x,2],[2,-y]])
minimize(2*x+y,[A<=0])

print x.sol, y.sol # sol: sqrt(2) , 2*sqrt(2)

print "eig(A) at optimum: ", np.linalg.eig(value(A))[0]

# Example 2: linear system
x0 = var(lb=0)
x1 = var(lb=0)

minimize(2*x0+x1*3,[x0+x1>=1])  

print x0.sol, x1.sol # 1 0


