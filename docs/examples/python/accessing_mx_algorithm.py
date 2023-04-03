#
#     MIT No Attribution
#
#     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
    assert v.is_dense()
    work[o[0]] = np.array(v)
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
      disp_in = ["work[" + str(a) + "]" for a in i]
      debug_str = print_operator(f.instruction_MX(k),disp_in)
      raise Exception('Unknown operation: ' + str(op) + ' -- ' + debug_str)

print('------')
print('Evaluated ' + str(f))
print('Expected: ', f.call(input_val))
print('Got:      ', output_val)
