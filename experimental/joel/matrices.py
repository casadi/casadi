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

# Create some symbolic variables
a = SX("a")
b = SX("b")
c = SX("c")
d = SX("d")

# Make a numpy matrix with numeric and symbolic variables
A = array([[a,0],[44,c]])
print "A = ", A

# Perform some operations (elementwise addition)
print "A + a = ", A + a

# Perform some operations (elementwise multiplication)
print "A * A = ", A * A

# Perform some operations (matrix multiplication)
print "dot(A,A) = ", dot(A,A)

# Create an SX sparse symbolic matrix (should be replaced by a typemap!!!)
B = SX(size(A,0),size(A,1))
for i in range(size(A,0)):
  for j in range(size(A,1)):
    B[i,j] = A[i,j]

# Print the SX
print "B = ", B

# Create a vector valued function
f = [sin(a)-3, exp(d-a)+a*a]
ffcn = SXFunction([[a,b,c,d]],[f]) # Note: single (vector valued input), single (vector valued) output

# Calculate the jacobian by source code transformation
jac_f = ffcn.jac()
print "jac_f = ", jac_f

# Number of non-zeros in the jacobian:
print "nnz(jac_f) = ", jac_f.size()

# Create a numpy array for continued symbolic processing (should be done with a typemap!)
J = zeros((jac_f.size1(),jac_f.size2()),dtype=SX)
for i in range(size(J,0)):
  for j in range(size(J,1)):
    J[i,j] = jac_f[i,j]
print "J = ", J





