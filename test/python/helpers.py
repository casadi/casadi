from casadi import *
from numpy import *
import unittest
import sys

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
  def message(self,s):
      print s
      sys.stdout.flush()
  def checkarray(self,zr,zt,name):
      """
      Checks for equality of two numpy matrices.
      The check uses dense form.
      
      zr - reference
      zt - to test
      
      name - a descriptor that will be included in error messages
      
      """
      if not(hasattr(zr,'shape')):
        zr=array([[zr]])
      if len(zr.shape)==1:
        zr=array(zr,ndmin=2).T
      self.assertEqual(zt.shape[0],zr.shape[0],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      self.assertEqual(len(zt.shape),len(zr.shape),"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      self.assertEqual(zt.shape[1],zr.shape[1],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
      for i in range(zr.shape[0]):
        for j in range(zr.shape[1]):
          self.assertAlmostEqual(zt[i,j],zr[i,j],10,"%s evaluation error. %s <-> %s" % (name, str(zt),str(zr)))


  def evaluationCheck(self,yt,yr,x,x0,name=""):
    """ General unit test for checking casadi evaluation against a reference solution.
    
        Checks if yr is the same as S/MXFunction(x,yt) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        yt: the test expression.
        yr: the reference solution: a numpy matrix.
        
        name - a descriptor that will be included in error messages
    """
    if (type(x)==list):
      sample = x[0]
    else :
      sample = x
      

    if isinstance(sample,SX) or isinstance(sample,SXMatrix):
      f = SXFunction(x,yt)
    else:
      f = MXFunction(x,yt)
      
    f.init()
    if not(type(x0)==list):
      x0=[x0]
    for i in range(len(x0)):
      try:
        f.setInput(x0[i],i)
      except:
         raise Exception("ERROR! Tried to set input with ", x0[i], " which is of type ", type(x0[i]))
    f.evaluate()
    zt = f.output(0).toArray()
    self.checkarray(yr,zt,name)
    
  def numpyEvaluationCheck(self,ft,fr,x,x0,name=""):
    """ General unit test for checking casadi versus numpy evaluation.
    
        Checks if 'fr(x0)' yields the same as S/MXFunction(x,[ft(x)]) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        ft: the test function. This function should operate on the casadi matrix x and return MX or SXMatrix.
        fr: the reference function. This function works on the numpy array x0.
        
        name - a descriptor that will be included in error messages
    """
    self.evaluationCheck([ft(x)],fr(x0),x,x0,name)
  
  def numpyEvaluationCheckPool(self,pool,x,x0,name=""):
    """ Performs a numpyEvaluationCheck for all members of a function pool"""
    for i in range(len(pool.numpyoperators)):
      self.numpyEvaluationCheck(pool.casadioperators[i],pool.numpyoperators[i],x,x0,name="\n I tried to apply %s (%s) from test case '%s' to numerical value %s. But the result returned: " % (str(pool.casadioperators[i]),pool.names[i],name, str(x0)))
