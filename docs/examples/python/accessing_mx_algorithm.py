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
Demonstration on how the algorithm of an MX function can be accessed and its operations can be transversed.
"""

from casadi import *
import numpy

# Create a function
a = MX.sym('a')
b = MX.sym('b',2)
c = MX.sym('c',2,2)
f = Function("f", [a,b,c], [3*mtimes(c,b)*a + b], ['a', 'b', 'c'], ['r'])

# Input values of the same dimensions as the above
input_val = [numpy.array([2.0]),\
             numpy.array([3.0,4.0]),\
             numpy.array([[5.0,1.0],[8.0,4.0]])]

# Output values to be calculated of the same dimensions as the above
output_val = [numpy.zeros(2)]

# Work vector
work = [ None for i in range(f.sz_w())]

# Loop over the algorithm
for k in range(f.n_instructions()):

  # Get the atomic operation
  op = f.instruction_id(k)

  o = f.instruction_output(k)
  i = f.instruction_input(k)

  if(op==OP_CONST):
    v = f.instruction_MX(k).to_DM()
    work[o[0]] = v
    print('work[{o[0]}] = {v}'.format(o=o,v=v))
  else:
    if op==OP_INPUT:
      work[o[0]] = input_val[i[0]]
      print('work[{o[0]}] = input[{i[0]}]            ---> {v}'.format(o=o,i=i,v=work[o[0]]))
    elif op==OP_OUTPUT:
      output_val[o[0]] = work[i[0]]
      print('output[{o[0]}] = work[{i[0]}]             ---> {v}'.format(o=o,i=i,v=output_val[o[0]]))
    elif op==OP_ADD:
      work[o[0]] = work[i[0]] + work[i[1]]
      print('work[{o[0]}] = work[{i[0]}] + work[{i[1]}]      ---> {v}'.format(o=o,i=i,v=work[o[0]]))
    elif op==OP_MUL:
      work[o[0]] = work[i[0]] * work[i[1]]
      print('work[{o[0]}] = work[{i[0]}] * work[{i[1]}]        ---> {v}'.format(o=o,i=i,v=work[o[0]]))
    elif op==OP_MTIMES:
      work[o[0]] = np.dot(work[i[1]], work[i[2]])+work[i[0]]
      print('work[{o[0]}] = work[{i[1]}] @ work[{i[2]}] + work[{i[0]}]        ---> {v}'.format(o=o,i=i,v=work[o[0]]))
    else:
      print('Unknown operation: ', op)

print('------')
print('Evaluated ' + str(f))
print('Expected: ', f.call(input_val))
print('Got:      ', output_val)
