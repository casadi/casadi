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
from time import time
from numpy import *
import casadi as c
import numpy as n


N = 5000

# We work with a sparse diaginal matrix with 5000 elements
x = DMatrix(sp_diag(N),5)

x_s = x.toCsr_matrix()

t = time()
dummy = n.dot(x_s.T,x_s)
print "Scipy pure = %.4f s" % (time()-t)

t = time()
dummy = c.mul(x.T,x)
print "CasADi pure = %.4f s" % (time()-t)

X = MX("X",x.sparsity())

t = time()
f = MXFunction([X],[c.mul(X.T,X)])
f.init()
print "CasADi MX wrapped init overhead = %.4f s" % (time()-t)
f.setInput(x)

t = time()
f.evaluate()
print "CasADi MX wrapped = %.4f s" % (time()-t)

t = time()
# Create the sparsity pattern for the matrix-matrix product
spres = x.sparsity().patternProduct(x.sparsity())
print "CasADi generating procuct sparsity pattern = %.4f s" % (time()-t)

t = time()
# Create the return object with correct sparsity
ret = c.DMatrix(spres)
print "CasADi allocating resulting = %.4f s" % (time()-t)

t = time()
# Carry out the matrix product
c.DMatrix.mul_no_alloc(x,x,ret)
print "CasADi actual multiplication = %.4f s" % (time()-t)
