#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

# For debugging
instr = f.instructions_sx()

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
      disp_in = ["work[" + str(a) + "]" for a in i]
      debug_str = print_operator(instr[k],disp_in)
      raise Exception('Unknown operation: ' + str(op) + ' -- ' + debug_str)

print('------')
print('Evaluated ' + str(f))
print('Expected: ', f.call(input_val))
print('Got:      ', output_val)
