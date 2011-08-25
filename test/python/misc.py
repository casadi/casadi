from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
import scipy.special
from scipy.linalg import expm

class Misctests(casadiTestCase):
    
  def test_issue179A(self):
    self.message('Regression test #179 (A)')
    x = SX("x")
    f = SXFunction([x], [2 * x])
    f.init()
    y = f.eval([x])[0].data()
    print y
    
  def test_issue179B(self):
    self.message('Regression test #179 (B)')
    def calc_sparsity():
      x = casadi.SX("x")
      f = casadi.SXFunction([[x]], [[x ** 2]])
      f.init()
      return f.jacSparsity()
    
    def print_sparsity():
        sparsity = calc_sparsity()
        print(sparsity) # Segfault
        
    print_sparsity()
    
  def test_sanity(self):
    DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99])
    self.assertRaises(RuntimeError,lambda : DMatrix(4,4,[1,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[0,2,2,12],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[-10,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[8,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[-3,2,1],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1,2],[0,2,2,3],[0.738,0.39,0.99]))
    self.assertRaises(RuntimeError,lambda : DMatrix(3,4,[1,2,1],[0,2,0,3],[0.738,0.39,0.99]))
    
    
    
if __name__ == '__main__':
    unittest.main()
    
