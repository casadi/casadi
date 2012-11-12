#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
# 
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
# 
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# 
# 
from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

#CasadiOptions.setCatchErrorsPython(False)


class FXtests(casadiTestCase):

  def test_call_empty(self):
    x = ssym("x",2)
    fsx = SXFunction([x,[]],[x])
    x = msym("x",2)
    fmx1 = MXFunction([x,MX()],[x])
    fmx2 = MXFunction([x,[]],[x])
    
    for f in [fsx,fmx1,fmx2]:
      f.init()
      f.evaluate(1,1)

      X = msym("X",2)
      F = f.call([X,X])[0]
      g = MXFunction([X],[F])
      g.init()

      g.evaluate(1,1)
    
    x = ssym("x",2)
    fsx = SXFunction([x],[x,[]])
    x = msym("x",2)
    fmx1 = MXFunction([x],[x,MX()])
    fmx2 = MXFunction([x],[x,[]])
    
    for f in [fsx,fmx1,]:
      f.init()
      f.evaluate(1,1)

      X = msym("X",2)
      F = f.call([X])[0]
      g = MXFunction([X],[F])
      g.init()

      g.evaluate(1,1)
  
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
    #self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0)) # This is now o.k. syntax
    B = DMatrix(1,2,2)
    #self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0)) # This is now o.k. syntax
    
  def test_issue304(self):
    self.message("regression test for #304") # this code used to segfault
    x = SX("x")

    f = SXFunction([x],[x**2,x**3])
    f.init()

    X = [MX("X")]

    z=f.call(X)

    g = MXFunction(X,[z[0]])
    g.init()

    g.expand([x])


      
if __name__ == '__main__':
    unittest.main()

