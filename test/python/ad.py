from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class ADtests(casadiTestCase):

  def setUp(self):
    x=SX("x")
    y=SX("y")
    z=SX("z")
    
    out=SXMatrix(5,1)
    out[0,0]=x
    out[2,0]=x+2*y**2
    out[4,0]=x+2*y**3+3*z**4
    

    inp=SXMatrix(5,1)
    inp[0,0]=x
    inp[2,0]=y
    inp[4,0]=z
    
    self.inputs = {
      "dense column vector": [[x,y,z]],
      "sparse column vector": [inp]
    }
    
    self.outputs = {
      "dense column vector": [[x,x+2*y**2,x+2*y**3+3*z**4]],
      "sparse column vector": [out]
    }
  
  def test_1(self):
    self.message("fwd AD of dense column vector wrt dense column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[1,4*y,0],[1,6*y**2,12*z**3]])
      
    f=SXFunction(self.inputs['dense column vector'],self.outputs['dense column vector'])
    f.init()
    f.input().set(n)
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.fwdSeed().set(d)
      f.evaluate(1,0)
      self.checkarray(f.fwdSens(),dot(reference(n),d),"AD")
    
  def test_2(self):
    self.message("adj AD of dense column vector wrt dense column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[1,4*y,0],[1,6*y**2,12*z**3]])
      
    f=SXFunction(self.inputs['dense column vector'],self.outputs['dense column vector'])
    f.init()
    f.input().set(n)
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.adjSeed().set(d)
      f.evaluate(0,1)
      self.checkarray(f.adjSens(),dot(reference(n).T,d),"AD")

  def test_3(self):
    self.message("fwd AD of sparse column vector wrt dense column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[0,0,0],[1,4*y,0],[0,0,0],[1,6*y**2,12*z**3]])
      
    f=SXFunction(self.inputs['dense column vector'],self.outputs['sparse column vector'])
    f.init()
    f.input().set(n)
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.fwdSeed().set(d)
      f.evaluate(1,0)
      self.checkarray(f.fwdSens(),dot(reference(n),d),"AD")

  def test_4(self):
    self.message("adj AD of sparse column vector wrt dense column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[0,0,0],[1,4*y,0],[0,0,0],[1,6*y**2,12*z**3]])
      
    f=SXFunction(self.inputs['dense column vector'],self.outputs['sparse column vector'])
    f.init()
    f.input().set(n)
    self.assertEqual(f.adjSeed().shape,(5,1),"adjSeed shape")
    self.assertEqual(f.adjSeed().size(),3,"adjSeed shape")
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.adjSeed().set(d)
      f.evaluate(0,1)
      self.checkarray(f.adjSens(),dot(reference(n).T,array(f.adjSeed())),"AD")

  def test_5(self):
    self.message("fwd AD of dense column vector wrt sparse column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0,0,0],[1,0,4*y,0,0],[1,0,6*y**2,0,12*z**3]])
    
    f=SXFunction(self.inputs['sparse column vector'],self.outputs['dense column vector'])
    f.init()
    f.input().set(n)
    self.assertEqual(f.fwdSeed().shape,(5,1),"adjSeed shape")
    self.assertEqual(f.fwdSeed().size(),3,"adjSeed shape")
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.fwdSeed().set(d)
      f.evaluate(1,0)
      self.checkarray(f.fwdSens(),dot(reference(n),array(f.fwdSeed())),"AD")

  def test6(self):
    self.message("adj AD of dense column vector wrt sparse column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0,0,0],[1,0,4*y,0,0],[1,0,6*y**2,0,12*z**3]])
    
    f=SXFunction(self.inputs['sparse column vector'],self.outputs['dense column vector'])
    f.init()
    f.input().set(n)
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.adjSeed().set(d)
      f.evaluate(0,1)
      self.checkarray(f.adjSens(),dot(reference(n).T,array(f.adjSeed())),"AD") 

  def test7(self):
    self.message("fwd AD of sparse column vector wrt sparse column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0,0,0],[0,0,0,0,0],[1,0,4*y,0,0],[0,0,0,0,0],[1,0,6*y**2,0,12*z**3]])
    
    f=SXFunction(self.inputs['sparse column vector'],self.outputs['sparse column vector'])
    f.init()
    f.input().set(n)
    self.assertEqual(f.fwdSeed().shape,(5,1),"adjSeed shape")
    self.assertEqual(f.fwdSeed().size(),3,"adjSeed shape")
    
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.adjSeed().set(d)
      f.evaluate(1,0)
      self.checkarray(f.fwdSens(),dot(reference(n),array(f.fwdSeed())),"AD") 
      
  def test8(self):
    self.message("adj AD of sparse column vector wrt sparse column vector")
    n=array([1.2,2.3,7])
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0,0,0],[0,0,0,0,0],[1,0,4*y,0,0],[0,0,0,0,0],[1,0,6*y**2,0,12*z**3]])
    
    f=SXFunction(self.inputs['sparse column vector'],self.outputs['sparse column vector'])
    f.init()
    f.input().set(n)
    self.assertEqual(f.adjSeed().shape,(5,1),"adjSeed shape")
    self.assertEqual(f.adjSeed().size(),3,"adjSeed shape")
    for d in [array([1,0,0]),array([0,2,0]),array([1.2,4.8,7.9])]:
      f.adjSeed().set(d)
      f.evaluate(0,1)
      self.checkarray(f.adjSens(),dot(reference(n).T,array(f.adjSeed())),"AD") 
           
if __name__ == '__main__':
    unittest.main()

