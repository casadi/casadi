from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class FXtests(casadiTestCase):
  
  def test_Parallelizer(self):
    self.message("Parallelizer")
    x = MX("x",2)
    y = MX("y")

    f = MXFunction([x,y],[sin(x) + y])
    f.init()
    
    #! Evaluate this function ten times in parallel
    p = Parallelizer([f]*2)
    
    for mode in ["openmp","serial"]:
      p.setOption("parallelization",mode)
      p.init()
      
      n1 = DMatrix([4,5])
      N1 = 3
      n2 = DMatrix([5,7])
      N2 = 8
      
      p.input(0).set(n1)
      p.input(1).set(N1)
      p.input(2).set(n2)
      p.input(3).set(N2)
      
      p.fwdSeed(0).setAll(0)
      p.fwdSeed(1).set(1)
      p.fwdSeed(2).set([1,0])
      p.fwdSeed(3).setAll(0)
      
      p.adjSeed(0).set([1,0])
      p.adjSeed(1).set([0,1])
      
      p.evaluate(1,1)

      self.checkarray(sin(n1)+N1,p.output(0),"output")
      self.checkarray(sin(n2)+N2,p.output(1),"output")
      self.checkarray(array([1,1]),p.fwdSens(0),"fwdSens")
      self.checkarray(array([cos(n2[0]),0]),p.fwdSens(1),"fwdSens")
      
      self.checkarray(array([cos(n1[0]),0]),p.adjSens(0),"adjSens")
      self.checkarray(1,p.adjSens(1),"adjSens")
      self.checkarray(array([0,cos(n2[1])]),p.adjSens(2),"adjSens")
      self.checkarray(1,p.adjSens(3),"adjSens")

  def test_MXFunctionSeed(self):
    self.message("MXFunctionSeed")
    x1 = MX("x",2)
    y1 = MX("y")
    x2 = MX("x",2)
    y2 = MX("y")
    p= MXFunction([x1,y1,x2,y2],[sin(x1) + y1,sin(x2) + y2])
    p.init()
    
    n1 = DMatrix([4,5])
    N1 = 3
    n2 = DMatrix([5,7])
    N2 = 8
    
    p.input(0).set(n1)
    p.input(1).set(N1)
    p.input(2).set(n2)
    p.input(3).set(N2)
    
    p.fwdSeed(0).setAll(0)
    p.fwdSeed(1).set(1)
    p.fwdSeed(2).set([1,0])
    p.fwdSeed(3).setAll(0)
    
    p.adjSeed(0).set([1,0])
    p.adjSeed(1).set([0,1])
    
    p.evaluate(1,1)

    self.checkarray(sin(n1)+N1,p.output(0),"output")
    self.checkarray(sin(n2)+N2,p.output(1),"output")
    self.checkarray(array([1,1]),p.fwdSens(0),"fwdSens")
    self.checkarray(array([cos(n2[0]),0]),p.fwdSens(1),"fwdSens")
    
    self.checkarray(array([cos(n1[0]),0]),p.adjSens(0),"adjSens")
    self.checkarray(1,p.adjSens(1),"adjSens")
    self.checkarray(array([0,cos(n2[1])]),p.adjSens(2),"adjSens")
    self.checkarray(1,p.adjSens(3),"adjSens")
    
  def test_Parallelizer2(self):
    self.message("Parallelizer")
    x = MX("x",2)
    y = MX("y")

    f = MXFunction([x,y],[sin(x) + y])
    f.init()
    
    #! Evaluate this function ten times in parallel
    pp = Parallelizer([f]*2)
    for mode in ["serial","openmp"]:
      pp.setOption("parallelization",mode)
      pp.init()
      
      x1 = MX("x",2)
      y1 = MX("y")
      x2 = MX("x",2)
      y2 = MX("y")
      p = MXFunction([x1,y1,x2,y2], pp.call([x1,y1,x2,y2]) )
      p.init()
      
      n1 = DMatrix([4,5])
      N1 = 3
      n2 = DMatrix([5,7])
      N2 = 8
      
      p.input(0).set(n1)
      p.input(1).set(N1)
      p.input(2).set(n2)
      p.input(3).set(N2)
      
      p.fwdSeed(0).setAll(0)
      p.fwdSeed(1).set(1)
      p.fwdSeed(2).set([1,0])
      p.fwdSeed(3).setAll(0)
      
      p.adjSeed(0).set([1,0])
      p.adjSeed(1).set([0,1])
      
      for i in range(4):
        p.adjSens(i).setAll(0)

      
      p.evaluate(1,1)

      self.checkarray(sin(n1)+N1,p.output(0),"output")
      self.checkarray(sin(n2)+N2,p.output(1),"output")
      self.checkarray(array([1,1]),p.fwdSens(0),"fwdSens")
      self.checkarray(array([cos(n2[0]),0]),p.fwdSens(1),"fwdSens")
      
      self.checkarray(array([cos(n1[0]),0]),p.adjSens(0),"adjSens")
      self.checkarray(1,p.adjSens(1),"adjSens")
      self.checkarray(array([0,cos(n2[1])]),p.adjSens(2),"adjSens")
      self.checkarray(1,p.adjSens(3),"adjSens")
    
  def test_ParallelizerMXCall(self):
    self.message("MX parallel call")
    x = MX("x",2)
    y = MX("y")

    f = MXFunction([x,y],[sin(x) + y])
    f.init()

    #! Evaluate this function ten times in parallel
    x1 = MX("x",2)
    y1 = MX("y")
    x2 = MX("x",2)
    y2 = MX("y")
    [[F1],[F2]] = f.call([[x1,y1],[x2,y2]])
    p = MXFunction([x1,y1,x2,y2],[F1,F2])
    p.init()
    
    n1 = DMatrix([4,5])
    N1 = 3
    n2 = DMatrix([5,7])
    N2 = 8
    
    p.input(0).set(n1)
    p.input(1).set(N1)
    p.input(2).set(n2)
    p.input(3).set(N2)
    
    p.fwdSeed(0).setAll(0)
    p.fwdSeed(1).set(1)
    p.fwdSeed(2).set([1,0])
    p.fwdSeed(3).setAll(0)
    
    p.adjSeed(0).set([1,0])
    p.adjSeed(1).set([0,1])
    
    p.evaluate(1,1)

    self.checkarray(sin(n1)+N1,p.output(0),"output")
    self.checkarray(sin(n2)+N2,p.output(1),"output")
    self.checkarray(array([1,1]),p.fwdSens(0),"fwdSens")
    self.checkarray(array([cos(n2[0]),0]),p.fwdSens(1),"fwdSens")


    self.checkarray(array([cos(n1[0]),0]),p.adjSens(0),"adjSens")
    self.checkarray(1,p.adjSens(1),"adjSens")
    self.checkarray(array([0,cos(n2[1])]),p.adjSens(2),"adjSens")
    self.checkarray(1,p.adjSens(3),"adjSens")

  def test_set_wrong(self):
    self.message("setter, wrong sparsity")
    x = SX("x")

    f = SXFunction([x],[x])
    f.init()
    A = DMatrix(1,1,4)
    f.getFwdSeed(A,0)
    A = DMatrix(1,1)
    self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0))
    B = DMatrix(1,2,2)
    self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0))
    
      
if __name__ == '__main__':
    unittest.main()

