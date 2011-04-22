from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class Matrixtests(casadiTestCase):
  def test_constructorlol(self):
    self.message("List of list constructor")
    a=DMatrix(array([[1,2,3],[4,5,6],[7,8,9]]))
    b=DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(a,b,"List of list constructor")
    
  def test_sum(self):
    self.message("sum")
    D=DMatrix([[1,2,3],[4,5,6],[7,8,9]])
    self.checkarray(c.sum(D),array([[12,15,18]]),'sum()')
    self.checkarray(c.sum(D,0),array([[12,15,18]]),'sum()')
    self.checkarray(c.sum(D,1),array([[6,15,24]]).T,'sum()')
    

    
if __name__ == '__main__':
    unittest.main()

