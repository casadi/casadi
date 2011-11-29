from casadi import *
from numpy import *
import unittest
import sys
from math import isnan, isinf

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
  def checkarray(self,zr,zt,name,failmessage=""):
      """
      Checks for equality of two numpy matrices.
      The check uses dense form.
      
      zr - reference
      zt - to test
      
      name - a descriptor that will be included in error messages
      
      """
      if isinstance(zr,tuple):
        zr=array([list(zr)])
      if isinstance(zt,tuple):
        zt=array([list(zt)])
      if not(hasattr(zt,'shape')) or len(zt.shape)==0:
        zt=array([[zt]])
      if not(hasattr(zr,'shape')) or len(zr.shape)==0:
        zr=array([[zr]])
      if len(zr.shape)==1:
        zr=array(zr,ndmin=2).T
      if len(zt.shape)==1:
        zt=array(zt,ndmin=2).T
        

      self.assertEqual(zt.shape[0],zr.shape[0],"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      self.assertEqual(len(zt.shape),len(zr.shape),"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      self.assertEqual(zt.shape[1],zr.shape[1],"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      for i in range(zr.shape[0]):
        for j in range(zr.shape[1]):
          if zt[i,j]==zr[i,j]:
            continue
          if (isnan(zt[i,j]) or isinf(zt[i,j])) and  (isinf(zt[i,j]) or isnan(zt[i,j])):
            continue
          self.assertAlmostEqual(zt[i,j],zr[i,j],10,"In %s: %s evaluation error.\n %s <->\n %s" % (name,failmessage, str(zt),str(zr)))

  def evaluationCheck(self,yt,yr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi evaluation against a reference solution.
    
        Checks if yr is the same as S/MXFunction(x,yt) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        yt: the test expression.
        yr: the reference solution: a numpy matrix.
        
        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    if (type(x)==list):
      sample = x[0]
    else :
      sample = x
      

    if isinstance(sample,SX) or isinstance(sample,SXMatrix):
      f = SXFunction(x,yt)
    else:
      f = MXFunction(x,yt)
    
    f.init()
    if (not(fmod is None)):
      f=fmod(f,x)
    if not(type(x0)==list):
      x0=[x0]
    
    if setx0 is None:
      setx0=x0
      
    if not(type(setx0)==list):
      setx0=[setx0]
      
    for i in range(len(x0)):
      try:
        f.input(i).set(setx0[i])
      except Exception as e:
         print f.input(i).shape
         raise e
         raise Exception("ERROR! Tried to set input with %s which is of type  %s \n%s" %(str(x0[i]), str(type(x0[i])),name))
    f.evaluate()
    zt = f.output(0).toArray()
    self.checkarray(yr,zt,name,failmessage)
    
  def numpyEvaluationCheck(self,ft,fr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi versus numpy evaluation.
    
        Checks if 'fr(x0)' yields the same as S/MXFunction(x,[ft(x)]) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        ft: the test function. This function should operate on the casadi matrix x and return MX or SXMatrix.
        fr: the reference function. This function works on the numpy array x0.
        
        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    fx=None
    frx=None
    try:
      fx=ft(x)
      frx=fr(x0)
    except Exception as e:
      print "Error calling functions in %s" % name
      raise e
    self.evaluationCheck([fx],frx,x,x0,name,failmessage,fmod=fmod,setx0=setx0)

  
  def numpyEvaluationCheckPool(self,pool,x,x0,name="",fmod=None,setx0=None):
    """ Performs a numpyEvaluationCheck for all members of a function pool"""
    for i in range(len(pool.numpyoperators)):
      self.numpyEvaluationCheck(pool.casadioperators[i],pool.numpyoperators[i],x,x0,"%s:%s" % (name,pool.names[i]),"\n I tried to apply %s (%s) from test case '%s' to numerical value %s. But the result returned: " % (str(pool.casadioperators[i]),pool.names[i],name, str(x0)),fmod=fmod,setx0=setx0)

