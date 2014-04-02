from casadi import *
import numpy

from helpers import *

class Symbolictests(casadiTestCase):

  @slow()
  @memory_heavy()
  def test_in_place_simplification(self):
    print "In place simplification"

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
        
    for a_s, b_s, xfunction in [(SX.sym("a"),SX.sym("b"),SXFunction),(MX.sym("a"),MX.sym("b"),MXFunction)]:
      print  xfunction
      i=0
      for op1 in operations_node():
        for op2 in operations_node():
           for op3 in operations_toplevel(op1,op2):
             i+= 1
             e = eval(op3,{"a": a_s, "b": b_s})
             f = xfunction([a_s,b_s],[e])
             #print i, xfunction, op1, op2, op3, " -> ", e
             f.init()
             f.setInput(numbers[0],0)
             f.setInput(numbers[1],1)
             f.evaluate()
             r = f.getOutput()
             try:
               num = eval(op3,{"a": numbers[0], "b":  numbers[1]})
             except ZeroDivisionError:
               num = None
             if num is None:
               self.assertTrue(numpy.isnan(r.at(0)) or numpy.isinf(r.at(0)))
               continue
             if (abs(num-r)>1e-10):
               print i , op3, " -> ", e
               self.assertTrue(False)

if __name__ == '__main__':
    unittest.main()
