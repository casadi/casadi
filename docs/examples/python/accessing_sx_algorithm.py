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
Demonstration on how the algorithm of an SX function can be accessed and its operations can be transversed.
"""

from casadi import *
import numpy

# Create a function
a = SX.sym('a')
b = SX.sym('b',2)
f = Function("f", [a,b], [2*a + b], ['a', 'b'], ['r'])

# Input values of the same dimensions as the above
input_val = [numpy.array([2.0]),\
             numpy.array([3.0,4.0])]

# Output values to be calculated of the same dimensions as the above
output_val = [numpy.zeros(2)]

# Work vector
work = numpy.zeros(f.sz_w())

# Loop over the algorithm
for k in range(f.n_instructions()):

  # Get the atomic operation
  op = f.instruction_id(k)
  o = f.instruction_output(k)
  i = f.instruction_input(k)

  if(op==OP_CONST):
    work[o[0]] = f.instruction_constant(k)
    print('work[', o[0], '] = ', f.instruction_constant(k))
  else:
    if op==OP_INPUT:
      work[o[0]] = input_val[i[0]][i[1]]
      print('work[', o[0], '] = input[', i[0], '][', i[1],  ']', '            ---> ' , work[o[0]])
    elif op==OP_OUTPUT:
      output_val[o[0]][o[1]] = work[i[0]]
      print('output[', o[0], '][', o[1], '] = work[', i[0], ']','             ---> ', output_val[o[0]][o[1]])
    elif op==OP_ADD:
      work[o[0]] = work[i[0]] + work[i[1]]
      print('work[', o[0], '] = work[', i[0], '] + work[', i[1], ']','        ---> ', work[o[0]])
    elif op==OP_MUL:
      work[o[0]] = work[i[0]] * work[i[1]]
      print('work[', o[0], '] = work[', i[0], '] * work[', i[1], ']','        ---> ', work[o[0]])
    else:
      print('Unknown operation: ', op)

print('------')
print('Evaluated ' + str(f))
print('Expected: ', f.call(input_val))
print('Got:      ', output_val)
