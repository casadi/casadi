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
import numpy

from helpers import *

class Symbolictests(casadiTestCase):

  @slow()
  @memory_heavy()
  def test_in_place_simplification(self):
    print("In place simplification")

    variables = ["a","b"]
    numbers = [1.3,0.7001]

    binary_operations = [lambda x,y :"(%s+%s)" % (x,y),lambda x,y: "(%s-%s)" % (x,y),lambda x,y :"(%s/%s)" % (x,y),lambda x,y: "(%s*%s)" % (x,y)]
    unary_operations = [lambda x: "(-%s)" % x, lambda x: x]


    """

       b1          b2               b3       b4
             op1                        op2
                           op3
    """


    def nodes():
      for b in variables:
        for u in unary_operations:
          yield u(b)

    def operations(a,b):
      for op in binary_operations:
        for u in unary_operations:
          yield u(op(a,b))
      for u in unary_operations:
        yield u(a)
        yield u(b)

    def operations_toplevel(a,b):
      for op in binary_operations:
        yield op(a,b)

    def operations_node():
      for b1 in nodes():
        for b2 in nodes():
          for op in operations(b1,b2):
            yield op
      for b1 in nodes():
        yield b1

    for a_s, b_s in [(SX.sym("a"),SX.sym("b")),(MX.sym("a"),MX.sym("b"))]:
      i=0
      for op1 in operations_node():
        for op2 in operations_node():
           for op3 in operations_toplevel(op1,op2):
             i+= 1
             e = eval(op3,{"a": a_s, "b": b_s})
             f = Function("f", [a_s,b_s],[e])
             #print i, op1, op2, op3, " -> ", e
             r = f(numbers[0], numbers[1])
             try:
               num = eval(op3,{"a": numbers[0], "b":  numbers[1]})
             except ZeroDivisionError:
               num = None
             if num is None:
               self.assertTrue(numpy.isnan(float(r)) or numpy.isinf(float(r)))
               continue
             if (abs(num-r)>1e-10):
               print(i , op3, " -> ", e)
               self.assertTrue(False)

if __name__ == '__main__':
    unittest.main()
