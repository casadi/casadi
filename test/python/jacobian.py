from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class Jacobiantests(casadiTestCase):
  def test_1(self):
    self.message("Jacobian of dense column vector")
    x=SX("x")
    y=SX("y")
    z=SX("z")
    n=array([1.2,2.3,7])
    f=SXFunction([[x,y,z]],[[x,x+2*y**2,x+2*y**3+3*z**4]])
    f.init()
    J=Jacobian(f,0,0)
    J.init()
    J.input().set(n)
    J.evaluate()
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[1,4*y,0],[1,6*y**2,12*z**3]])
    
    self.checkarray(reference(n),array(J.output()),"Jacobian of dense column vector")
   
  def test_2(self):
    return #bug
    self.message("Jacobian of sparse column vector")
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
    J=Jacobian(f,0,0)
    J.init()
    J.input().set(n)
    J.evaluate()
    def reference(n):
      x=n[0]
      y=n[1]
      z=n[2]
      return array([[1,0,0],[0,0,0],[1,4*y,0],[0,0,0],[1,6*y**2,12*z**3]])
    
    self.checkarray(reference(n),array(J.output()),"Jacobian of dense column vector")
    
if __name__ == '__main__':
    unittest.main()

