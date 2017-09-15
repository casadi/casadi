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
# -*- coding: utf-8 -*-
from casadi import *

# Tests the finite difference implementation in CasADi
# Joel Andersson, UW Madison 2017

# The following function is defined for the interval [0,1] and gives NaN outside it
x = SX.sym('x');
r = if_else(x<=1, if_else(x>=0, 2*x+cos(x), log(x)), log(1-x))

# Test points in the domain, on the boundary and outside it
x_test = [-1., 0., 0.4, 1., 1.4]
print("%20s" % "x0", end='')
for x0 in x_test: print("%15g " % x0, end='')
print()

# Try all supported schemes
for fd_method in ['forward','backward','central','smoothing']:
    print("%20s" % fd_method, end='')
    # Construct function and differentiate
    opts = dict()
    opts["enable_fd"]=True # enable finite differencing
    opts["enable_forward"]=False # disable forward mode AD
    opts["enable_reverse"]=False # disable reverse mode AD
    opts["enable_jacobian"]=False # disable AD by calculating full Jacobian
    opts["fd_method"]=fd_method # specify FD scheme
    f = Function('f', [x], [r], ['x'], ['r'], opts)
    fwd_f = f.forward(1)
    # Evaluate
    for x0 in x_test:
        r0 = f(x0)
        d = fwd_f(x0, r0, 1.)
        print("%15g " % float(d), end='')
    print()
