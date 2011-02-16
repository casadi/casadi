from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class Jacobiantests(casadiTestCase):
  def test_1(self):
    self.message("Jacobian of dense column vector wrt dense column vector")
    x=SX("x")
    y=SX("y")
    z=SX("z")
    n=array([1.2,2.3,7])
    f=SXFunction([[x,y,z]],[[x,x+2*y**2,x+2*y**3+3*z**4]])
    f.init()
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[1,4*y,0],[1,6*y**2,12*z**3]])
    for mode in ['forward','adjoint']:
      J=Jacobian(f,0,0)
      J.setOption("ad_mode",mode)
      J.init()
      J.input().set(n)
      J.evaluate()
      self.checkarray(reference(n),array(J.output()),"Jacobian")
   
  def test_2(self):
    return #bug
    self.message("Jacobian of sparse column vector wrt dense column vector")
    x=SX("x")
    y=SX("y")
    z=SX("z")
    n=array([1.2,2.3,7])
    out=SXMatrix(5,1)
    out[0,0]=x
    out[2,0]=x+2*y**2
    out[4,0]=x+2*y**3+3*z**4
    
    f=SXFunction([[x,y,z]],[out])
    f.init()
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[0,0,0],[1,4*y,0],[0,0,0],[1,6*y**2,12*z**3]])
    
    for mode in ['forward','adjoint']:
      J=Jacobian(f,0,0)
      J.setOption("ad_mode",mode)
      J.init()
      J.input().set(n)
      J.evaluate()
      self.checkarray(reference(n),array(J.output()),"Jacobian")

  def test_3(self):
    return #glibc bug
    self.message("Jacobian of dense column vector wrt sparse column vector")
    x=SX("x")
    y=SX("y")
    z=SX("z")
    n=array([1.2,2.3,7])
    inp=SXMatrix(5,1)
    inp[0,0]=x
    inp[2,0]=y
    inp[4,0]=z
    
    f=SXFunction([inp],[[x,x+2*y**2,x+2*y**3+3*z**4]])
    f.init()
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0,0,0],[1,0,4*y,0,0],[1,0,6*y**2,0,12*z**3]])
    for mode in ['forward','adjoint']:
      J=Jacobian(f,0,0)
      J.setOption("ad_mode",mode)
      J.init()
      J.input().set(n)
      J.evaluate()
      self.checkarray(reference(n),array(J.output()),"Jacobian")
      
  def test_bugshape(self):
    return
    x=SX("x")
    y=SX("y")

    inp=SXMatrix(5,1)
    inp[0,0]=x
    inp[3,0]=y

    f=SXFunction([inp],[[x+y,x,y]])
    f.init()
    J=Jacobian(f,0,0)
    J.setOption("ad_mode","forward")
    J.init()
    J.input().set([2,7])
    J.evaluate()

    self.assertEqual(f.output().size1(),3,"Jacobian shape bug")
    self.assertEqual(f.output().size2(),1,"Jacobian shape bug")

    
  def test_bugglibc(self):
    return
    self.message("Code that used to throw a glibc error")
    x=SX("x")
    y=SX("y")

    inp=SXMatrix(5,1)
    inp[0,0]=x
    inp[3,0]=y

    f=SXFunction([inp],[[x+y,x,y]])
    f.init()
    J=Jacobian(f,0,0)
    J.setOption("ad_mode","forward")
    J.init()
    J.input().set([2,7])
    J.evaluate()

    f=SXFunction([inp],[[x+y,x,y]])
    f.init()
    print f.input().shape
    J=Jacobian(f,0,0)
    
if __name__ == '__main__':
    unittest.main()

