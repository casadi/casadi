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
from optoy import *

# Simple unconstrained rosenbrock problem
x = var()
y = var()

cost = minimize((1-x)**2+100*(y-x**2)**2)
print "cost = ", cost
print "sol = ", x.sol, y.sol

# Bounds
x.lb = 2
y = var(ub=6)

print "cost = ", minimize((1-x)**2+100*(y-x**2)**2)
print "sol = ", x.sol, y.sol

# Constraints
x = var()
y = var()

print "cost = ", minimize((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])
print "sol = ", x.sol, y.sol

# Matrix symbols
xy = var(2)

print "cost = ", minimize((1-xy[0])**2+100*(xy[1]-xy[0]**2)**2)
print "sol = ", xy.sol

# Adding a parameter
p = par()
p.value = 100

print "cost = ", minimize((1-x)**2+p*(y-x**2)**2)
print "sol = ", x.sol, y.sol

# Subexpressions
x = var()
y = var()

e = (y-x**2)**2

cost = minimize((1-x)**2+100*e)
print "cost = ", cost
print "sol = ", x.sol, y.sol
print "e = ", value(e)
