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
"""
Demonstration on how the algorithm of an SXFunction can be accessed and its operations can be transversed.
"""

from casadi import *
import numpy

a = SX.sym('a')
b = SX.sym('b',2)

# Input expressions
input_ex = [a,b]

# Input values of the same dimensions as the above
input_val = [numpy.array([2.0]),numpy.array([3.0,4.0])]

# Output expressions
output_ex = [2*a + b]

# Output values to be calculated of the same dimensions as the above
output_val = [numpy.zeros(2)]

# Create a function
f = SXFunction(input_ex,output_ex)
f.init()

# Work vector
work = numpy.zeros(f.getWorkSize())

# Loop over the algorithm
for i in range(f.getAlgorithmSize()):
  
  # Get the atomic operation
  op = f.getAtomicOperation(i)
  
  if(op==OP_CONST):
    work[f.getAtomicOutput(i)] = f.getAtomicInputReal(i)
    print 'work[', f.getAtomicOutput(i), '] = ', f.getAtomicInputReal(i)
  else:
    i1 = f.getAtomicOutput(i)
    i2,i3 = f.getAtomicInput(i)
    if op==OP_INPUT:
      work[i1] = input_val[i2][i3]
      print 'work[', i1, '] = input[', i2, '][', i3,  ']', '                ---> ' , work[i1]
    elif op==OP_OUTPUT:
      output_val[i1][i3] = work[i2]
      print 'output[', i1, '][', i3, '] = work[', i2, ']','             ---> ', output_val[i1][i3]
    elif op==OP_ADD:
      work[i1] = work[i2] + work[i3]
      print 'work[', i1, '] = work[', i2, '] + work[', i3, ']','        ---> ', work[i1]
    elif op==OP_MUL:
      work[i1] = work[i2] * work[i3]
      print 'work[', i1, '] = work[', i2, '] * work[', i3, ']','        ---> ', work[i1]
    else:
      print 'Unknown operation: ', op

print '------'
print "Evaluated function: "
print output_ex, ' = ', output_val
print 'where ', input_ex, ' = ', input_val
      
