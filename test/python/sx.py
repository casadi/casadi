from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class SXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: sqrt(x[0]),sqrt,"sqrt")
    self.pool.append(lambda x: sin(x[0]),sin,"sin")
    self.pool.append(lambda x: cos(x[0]),cos,"cos")
    self.pool.append(lambda x: tan(x[0]),tan,"tan")
    self.pool.append(lambda x: arctan(x[0]),arctan,"arctan")
    self.pool.append(lambda x: arcsin(x[0]),arcsin,"arcsin")
    self.pool.append(lambda x: arccos(x[0]),arccos,"arccos")
    self.pool.append(lambda x: exp(x[0]),exp,"exp")
    self.pool.append(lambda x: log(x[0]),log,"log")
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0")
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2")
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: floor(x[0]),floor,"floor")
    self.pool.append(lambda x: ceil(x[0]),ceil,"ceil")
    #self.pool.append(lambda x: erf(x[0]),erf,"erf") # numpy has no erf
    
    
  def test_scalarSX(self):
      x=SXMatrix("x")
      x0=0.738
      
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarSX")
      
  def test_SXMAtrix(self):
      x=SXMatrix("x",3,2)
      x0=array([[0.738,0.2],[ 0,0.39 ],[0.99,0.999999]])
      
      #self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarSX")
  
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

