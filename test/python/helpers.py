from casadi import *
from numpy import *
import unittest

class FunctionPool:
  def __init__(self):
    self.numpyoperators=[]
    self.casadioperators=[]
    self.names=[]
  def append(self,cas,num,name=""):
    self.casadioperators.append(cas)
    self.numpyoperators.append(num)
    self.names.append(name)

class casadiTestCase(unittest.TestCase):
  def checkarray(self,zr,zt,name):
      """
      Checks for equality of two numpy matrices.
      The check uses dense form.
      
      zr - reference
      zt - to test
      
      name - a descriptor that will be included in error messages
      
      """
      zr=array(zr,ndmin=2)
      self.assertEqual(zt.shape[0],zr.shape[0],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      self.assertEqual(len(zt.shape),len(zr.shape),"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      self.assertEqual(zt.shape[1],zr.shape[1],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      for i in range(zr.shape[0]):
        for j in range(zr.shape[1]):
          self.assertAlmostEqual(zt[i,j],zr[i,j],10,"%s evaluation error. %s <-> %s" % (name, str(zt),str(zr)))



  def numpyEvaluationCheck(self,ft,fr,x,x0,name=""):
    """ General unit test for checking casadi versus numpy evaluation.
    
        Checks if 'fr(x0)' yields the same as S/MXFunction(x,[y]) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        ft: the test function. This function should operate on the casadi matrix x and return MX or SXMatrix.
        fr: the reference function. This function works on the numpy array x0.
        
        name - a descriptor that will be included in error messages
    """
    if (type(x)==list):
      sample = x[0]
    else :
      sample = x
      
    y=ft(x)
    if isinstance(sample,SX) or isinstance(sample,SXMatrix):
      f = SXFunction(x,[y])
    else:
      f = MXFunction(x,[y])
      
    f.init()
    f.setInput(x0,0)
    f.evaluate()
    yt = f.output(0).getArray()
    yr = fr(x0)
    self.checkarray(yr,yt,name)
  
  def numpyEvaluationCheckPool(self,pool,x,x0,name=""):
    """ Performs a numpyEvaluationCheck for all members of a function pool"""
    for i in range(len(pool.numpyoperators)):
      self.numpyEvaluationCheck(pool.casadioperators[i],pool.numpyoperators[i],x,x0,name="\n I tried to apply %s (%s) from test case '%s' to numerical value %s. But the result returned: " % (str(pool.casadioperators[i]),pool.names[i],name, str(x0)))
