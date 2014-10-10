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
#! CasADi tutorial 3
#! ==================
#! This tutorial file explains the use of CasADi's MX in a python context.
from numpy import *
from casadi import *
#! Let's start with some algebra
x = MX.sym("x",2,3)
y = MX.sym("y",3,2)
print x
for i in range(6):
	print x.nz[i]

for i in range(2):
	for j in range(3):
		print "x[%d,%d] = %s" % (i,j,str(x[i,j]))
		
print x[1,1]
print x.nz[3] # Note that index is flattened. x[0,0] is illegal.
print norm_2(x)
z= mul(x,y)
print z
#! Note how the operations on MXes are lazy on the matrix level.
#! Any elementwise logic is postponed until evaluation demands it.
#! Just like, SXFunction, MXFunction can be single or multi input/output.
#! The only allowed input/output primitive is MX.
f = MXFunction([x,y],[z])


#! Evaluation
#! ----------------
f.init()
f.setInput([1,2,3,4,5,6],0);
f.setInput([1,3,0,6,0,9],1);
f.evaluate()

print f.getOutput()
#! Note how this result is related to a numpy approach:
a=matrix(f.getInput(0)).reshape(3,2)
b=matrix(f.getInput(1)).reshape(2,3)
print a.T*b.T
#! Jacobian
#! -------------
#! BUG error
#!print f.jacobian(0)
#!type(f.hessian())

#! Numerical matrices
#! ------------------------------
X = MX(DMatrix(array([[1,2,3],[4,5,6]])))
print X
print outer_prod(X,X)
print MX(DMatrix([1,2,3]).T)
print MX([1,2,3])
#! As before, evaluation is lazy on the matrix level
Y = MX.sym("Y")
f = MXFunction([Y],[X])
f.init()
f.setInput([2])
f.evaluate()
print f.getOutput()

#! Element assignement
#! -------------------
X = MX.sym("x",2,2)
X[0,0]=MX(5)
print X
 
