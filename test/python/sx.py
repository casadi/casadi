from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *

class SXtests(unittest.TestCase):

  def setUp(self):
    pass
    
  def test_SX1(self):
    fun=lambda x,y: [x+y,x*y,x**2+y**3]
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[fun(x,y)])
    f.init()
    L=[2,3]
    f.setInput(L)
    f.evaluate()
    z=f.output(0).getArray()
    zr=fun(*L)
    for i in range(3):
      self.assertAlmostEqual(z[i], zr[i],10,'SXfunction output in correct')
    
    J=f.jacobian()
    J.init()
    J.setInput(L)
    J.evaluate()
    J=J.output(0).getArray()
    Jr=matrix([[1,1],[3,2],[4,27]])
    for i in range(3):
        for j in range(2):
          self.assertAlmostEqual(J[i,j], Jr[i,j],10,'SXfunction jacobian in correct')
          
  def test_SX2(self):
    fun = lambda x,y: [3-sin(x*x)-y, sqrt(y)*x]
    # variables
    x = SX("x")
    y = SX("y")

    # Create function
    f = fun(x,y)
    self.assertEqual(str(f),'[((3-sin((x*x)))-y), (sqrt(y)*x)]','SX representation is wrong')
    fcn = SXFunction([[x,y]],[f])

    # Set some options
    fcn.setOption("name","f")
    fcn.setOption("ad_order",1)

    self.assertEqual(str(fcn),'sx function("f")','SX representation is wrong')
    # Initialize the function for numerical calculations
    fcn.init()

    # Pass inputs
    L=[2,3]
    fcn.setInput(L)

    # Pass forward seed
    sF=[1,0]
    fcn.setFwdSeed(sF);

    sA=[1,0]
    # Pass adjoint seed
    fcn.setAdjSeed(sA)

    # Evaluate numerically
    fcn.evaluate(1,1)

    # Get the results
    res = fcn.getOutput()
    self.assertAlmostEqual(res[0], fun(*L)[0],10,'SXfunction evaluation wrong')
    self.assertAlmostEqual(res[1], fun(*L)[1],10,'SXfunction evaluation wrong')

    fsens = fcn.getFwdSens()
    e=1e-8
    p=fun((L[0]+e*sF[0]),(L[1]+e*sF[1]))
    fsensn=[(p[0]-res[0])/e,(p[1]-res[1])/e]
    self.assertAlmostEqual(fsensn[0], fsens[0],6,'SXfunction forward mode evaluation wrong')
    self.assertAlmostEqual(fsensn[1], fsens[1],6,'SXfunction forward mode evaluation wrong')
    
    J=fcn.jacobian()
    J.init()
    J.setInput(L)
    J.evaluate()
    J=J.output(0).getArray()
    
    fsensJ=dot(J,array(sF))
    self.assertAlmostEqual(fsensJ[0], fsens[0],10,'SXfunction forward mode evaluation wrong')
    self.assertAlmostEqual(fsensJ[1], fsens[1],10,'SXfunction forward mode evaluation wrong')
    
    asens = fcn.getAdjSens()
    asensJ=dot(J.T,array(sA))
    self.assertAlmostEqual(asensJ[0], asens[0],10,'SXfunction adjoint mode evaluation wrong')
    self.assertAlmostEqual(asensJ[1], asens[1],10,'SXfunction adjoint mode evaluation wrong')

    # evaluate symbolically
    fsubst = fcn.eval([[1-y,1-x]])
    #print "fsubst = ", fsubst
      
if __name__ == '__main__':
    unittest.main()

