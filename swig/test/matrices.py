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

# Create an SXMatrix sparse symbolic matrix (should be replaced by a typemap!!!)
B = SXMatrix(size(A,0),size(A,1))
for i in range(size(A,0)):
  for j in range(size(A,1)):
    B[i,j] = A[i,j]

# Print the SXMatrix
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





