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
from numpy import *
import matplotlib.pyplot as plt

A = array( [[ 0.0,  0.0, 6.0, 7.0, 8.0, 0.0 ],
            [ 0.0,  0.0, 0.0, 9.0,10.0,11.0 ],
            [ 0.0,  1.0, 2.0, 0.0, 0.0, 0.0 ],
            [ 14.0, 0.0, 0.0, 0.0,12.0,13.0 ],
            [  0.0, 3.0, 4.0, 5.0, 0.0, 0.0 ],
            [ 16.0, 0.0, 0.0, 0.0, 0.0,15.0 ]])

b = array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

x = linalg.solve(A,b)
print "x  = ", x

cA = DMatrix(A)
makeSparse(cA)

cb = DMatrix(b)
cx = casadi.solve(cA,cb)
print "cx = ", array(trans(cx))[0]

A = array( [[ 0.0,  1.0 ],
            [ 1.0,  0.0 ]] )
b = array([1.0, 2.0])

x = linalg.solve(A,b)
print "x  = ", x

cA = DMatrix(A)
makeSparse(cA)

cb = DMatrix(b)
cx = casadi.solve(cA,cb)
print "cx = ", array(trans(cx))[0]
