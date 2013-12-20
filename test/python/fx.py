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
      
      p.setInput(n1,0)
      p.setInput(N1,1)
      p.setInput(n2,2)
      p.setInput(N2,3)
      
      p.setFwdSeed(0,0)
      p.setFwdSeed(1,1)
      p.setFwdSeed([1,0],2)
      p.setFwdSeed(0,3)
      
      p.setAdjSeed([1,0],0)
      p.setAdjSeed([0,1],1)
      
      p.evaluate(1,1)

      self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
      self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
      self.checkarray(array([1,1]),p.getFwdSens(0),"fwdSens")
      self.checkarray(array([cos(n2[0]),0]),p.getFwdSens(1),"fwdSens")
      
      self.checkarray(array([cos(n1[0]),0]),p.getAdjSens(0),"adjSens")
      self.checkarray(1,p.getAdjSens(1),"adjSens")
      self.checkarray(array([0,cos(n2[1])]),p.getAdjSens(2),"adjSens")
      self.checkarray(1,p.getAdjSens(3),"adjSens")

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
    
    p.setInput(n1,0)
    p.setInput(N1,1)
    p.setInput(n2,2)
    p.setInput(N2,3)
    
    p.setFwdSeed(0,0)
    p.setFwdSeed(1,1)
    p.setFwdSeed([1,0],2)
    p.setFwdSeed(0,3)
    
    p.setAdjSeed([1,0],0)
    p.setAdjSeed([0,1],1)
    
    p.evaluate(1,1)

    self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
    self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
    self.checkarray(array([1,1]),p.getFwdSens(0),"fwdSens")
    self.checkarray(array([cos(n2[0]),0]),p.getFwdSens(1),"fwdSens")
    
    self.checkarray(array([cos(n1[0]),0]),p.getAdjSens(0),"adjSens")
    self.checkarray(1,p.getAdjSens(1),"adjSens")
    self.checkarray(array([0,cos(n2[1])]),p.getAdjSens(2),"adjSens")
    self.checkarray(1,p.getAdjSens(3),"adjSens")
    
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
      
      p.setInput(n1,0)
      p.setInput(N1,1)
      p.setInput(n2,2)
      p.setInput(N2,3)
      
      p.setFwdSeed(0,0)
      p.setFwdSeed(1,1)
      p.setFwdSeed([1,0],2)
      p.setFwdSeed(0,3)
      
      p.setAdjSeed([1,0],0)
      p.setAdjSeed([0,1],1)
      
      for i in range(4):
        p.setAdjSens(0,i)

      
      p.evaluate(1,1)

      self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
      self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
      self.checkarray(array([1,1]),p.getFwdSens(0),"fwdSens")
      self.checkarray(array([cos(n2[0]),0]),p.getFwdSens(1),"fwdSens")
      
      self.checkarray(array([cos(n1[0]),0]),p.getAdjSens(0),"adjSens")
      self.checkarray(1,p.getAdjSens(1),"adjSens")
      self.checkarray(array([0,cos(n2[1])]),p.getAdjSens(2),"adjSens")
      self.checkarray(1,p.getAdjSens(3),"adjSens")
    
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
    
    p.setInput(n1,0)
    p.setInput(N1,1)
    p.setInput(n2,2)
    p.setInput(N2,3)
    
    p.setFwdSeed(0,0)
    p.setFwdSeed(1,1)
    p.setFwdSeed([1,0],2)
    p.setFwdSeed(0,3)
    
    p.setAdjSeed([1,0],0)
    p.setAdjSeed([0,1],1)
    
    p.evaluate(1,1)

    self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
    self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
    self.checkarray(array([1,1]),p.getFwdSens(0),"fwdSens")
    self.checkarray(array([cos(n2[0]),0]),p.getFwdSens(1),"fwdSens")


    self.checkarray(array([cos(n1[0]),0]),p.getAdjSens(0),"adjSens")
    self.checkarray(1,p.getAdjSens(1),"adjSens")
    self.checkarray(array([0,cos(n2[1])]),p.getAdjSens(2),"adjSens")
    self.checkarray(1,p.getAdjSens(3),"adjSens")

  def test_set_wrong(self):
    self.message("setter, wrong sparsity")
    x = SX("x")

    f = SXFunction([x],[x])
    f.init()
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
  
  def test_jacobian(self):
    x = ssym("x",3,1)
    y = ssym("y",2,1)

    f = SXFunction([x,y],[x**2,y,x*y[0]])
    f.init()

    g = f.jacobian(0,0)

    self.assertEqual(g.getNumInputs(),f.getNumInputs())
    self.assertEqual(g.getNumOutputs(),f.getNumOutputs()+1)

  def test_xfunction(self):
    x = ssym("x",3,1)
    y = ssym("y",2,1)
    
    f = SXFunction([x,y],[x**2,y,x*y[0]])
    f.init()
    
    f.setInput([0.1,0.7,1.3],0)
    f.setInput([7.1,2.9],1)
    
    X = msym("x",3,1)
    Y = msym("y",2,1)
    
    F = MXFunction([X,Y],[X**2,Y,X*Y[0]])
    F.init()
    
    F.setInput([0.1,0.7,1.3],0)
    F.setInput([7.1,2.9],1)
    
    self.checkfx(f,F,sens_der=False,evals=False)
    
  
  @memory_heavy()
  def test_jacobians(self):
  
    x = ssym("x")
    
    self.assertEqual(jacobian(5,x).size(),0)
    
    
    def test(sp):
      x = ssym("x",sp.size2())
      sp2 = jacobian(mul(DMatrix(sp,1),x),x).sparsity()
      self.checkarray(sp.col(),sp2.col());
      self.checkarray(sp.rowind(),sp2.rowind());   

    for i in range(5):
      test(sp_tril(i))
      test(sp_tril(i).T)
      test(sp_dense(i,i))
      test(sp_diag(i))
    
    for i in [63,64,65,127,128,129]:
      d = sp_diag(i)
      test(d)
      
      test(d + sp_rowcol([0],[5],i,i))
      
      b = sp_band(i,-1) + sp_band(i,1)
      test(b + sp_rowcol([0],[5],i,i))
      
    m = IMatrix(sp_diag(129),1)
    m[:50,0] = 1
    m[60:,0] = 1
    m[6:9,6] = 1
    m[9,9:12] = 1
    
    sp = m[:,:120].sparsity()
    
    test(sp)
    #test(sp.T)
    
    m = IMatrix(sp_diag(64),1)
    m[:50,0] = 1
    m[60:,0] = 1

    sp = m.T[:40,:].sparsity()
    test(sp)
    test(sp.T)
    
    sp = m[:40,:].sparsity()
    test(sp)
    test(sp.T)
    
    sp = m.T[:20,:].sparsity()
    test(sp)
    test(sp.T)

    sp = m[:20,:].sparsity()
    test(sp)
    test(sp.T)
    
    for i in [63,64,65,127,128,129]:
      test(sp_tril(i))
      test(sp_tril(i).T)
    
    for n in [63,64,65,127,128,129]:
      for m in [63,64,65,127,128,129]:
        print (n,m)
        sp = sp_dense(n,m)
        
        test(sp)
        
        random.seed(0)
        
        I = IMatrix(sp,1)
        for i in range(n):
          for j in range(m):
            if random.random()<0.5:
              I[i,j] = 0
        makeSparse(I)
        
        sp_holes = I.sparsity()
        
        test(sp_holes)
        
        z = IMatrix(sp_holes.shape[0],sp_holes.shape[1])
        
        R = 5
        v = []
        for r in range(R):
          h = [z]*5
          h[r] = I
          v.append(horzcat(h))
        d = vertcat(v)
        
        test(d.sparsity())
        
  @memory_heavy()
  def test_hessians(self):
    def test(sp):
      x = ssym("x",sp.size2())
      self.assertTrue(sp==sp.transpose())
      f = SXFunction([x],[mul([x.T,DMatrix(sp,1),x])])
      f.init()
      J = f.hessian()
      J.init()
      sp2 = J.output().sparsity()
      self.checkarray(sp.col(),sp2.col())
      self.checkarray(sp.rowind(),sp2.rowind())
      
    A = IMatrix([[1,1,0,0,0,0],[1,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    makeSparse(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    makeSparse(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    makeSparse(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    makeSparse(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,0,0,1,0,1]])
    makeSparse(A)
    C = A.sparsity()
      
    test(C)
    
    
    for i in [63,64,65,100,127,128,129]:
      d = sp_diag(i)
      test(d)
      
      test(d + sp_rowcol([0],[5],i,i) + sp_rowcol([5],[0],i,i))
      
      b = sp_band(i,-1) + sp_band(i,1)
      test(b)
      test(b + sp_rowcol([0],[5],i,i) + sp_rowcol([5],[0],i,i))
      
      d = sp_dense(i,i)
      test(d)
      
      d = sp_diag(i) + sp_triplet(i,i,[0]*i,range(i))+sp_triplet(i,i,range(i),[0]*i)
      test(d)


      sp = sp_dense(i,i)
        
      random.seed(0)
      
      I = IMatrix(sp,1)
      for ii in range(i):
        for jj in range(i):
          if random.random()<0.5:
            I[ii,jj] = 0
            I[jj,ii] = 0
      makeSparse(I)
      
      sp_holes = I.sparsity()
      
      test(sp_holes)
      
      z = IMatrix(sp_holes.shape[0],sp_holes.shape[1])
      
      R = 5
      v = []
      for r in range(R):
        h = [z]*5
        h[r] = I
        v.append(horzcat(h))
      d = vertcat(v)
      
      test(d.sparsity())
      
  def test_getOutput(self):
    x = ssym("x",2)
    
    f = SXFunction(daeIn(x=x),daeOut(ode=x))
    f.init()
    f.setInput([1,2])
    f.evaluate()
    a = f.getOutput()
    b = f.getOutput(0)
    c = f.getOutput("ode")
    ar = f.output()
    br = f.output(0)
    cr = f.output("ode")
    self.checkarray(a,DMatrix([1,2]))
    self.checkarray(b,DMatrix([1,2]))
    self.checkarray(c,DMatrix([1,2]))
    self.checkarray(ar,DMatrix([1,2]))
    self.checkarray(br,DMatrix([1,2]))
    self.checkarray(cr,DMatrix([1,2]))
    f.setInput([3,4])
    f.evaluate()
    self.checkarray(a,DMatrix([1,2]))
    self.checkarray(b,DMatrix([1,2]))
    self.checkarray(c,DMatrix([1,2]))
    self.checkarray(ar,DMatrix([3,4]))
    self.checkarray(br,DMatrix([3,4]))
    self.checkarray(cr,DMatrix([3,4]))
    
  def test_customIO(self):
    
    myOut = IOScheme(["foo","bar"])

    x = ssym("x")
    
    with self.assertRaises(Exception):
      myOut(baz=x)
    
    f = SXFunction([x],myOut(foo=x*x,bar=x))
    f.init()
    
    f.setInput(12,0)
    f.evaluate()
    self.checkarray(DMatrix([144]),f.output("foo"))
    self.checkarray(DMatrix([12]),f.output("bar"))

    
    with self.assertRaises(Exception):
      f.output("baz")
      
    ret = f.eval([12])
    self.checkarray(myOut(ret,"foo")[0],DMatrix([144]))
    self.checkarray(myOut(ret,"bar")[0],DMatrix([12]))
    with self.assertRaises(Exception):
      self.checkarray(myOut(ret,"baz")[0],DMatrix([12]))
     
      
  def test_unknown_options(self):
    x = ssym("x")
    f = SXFunction([x],[x])
    f.init()
    
    with self.assertRaises(Exception):
      f.setOption({"fooo": False},False)
    
    f.setOption({"fooo": False},True)
    
    f.setOption({"name": "abc"},False)
    self.assertTrue(f.getOption("name")=="abc")
    f.setOption({"name": "def"},True)
    self.assertTrue(f.getOption("name")=="def")
    
  def test_CustomFunctionDefault(self):
    x = msym("x")
    y = msym("y")

    @pyevaluate
    def fun(f,nfwd,nadj):
      # sin(x+3*y)
      
      print (nfwd,nadj)
      
      x = f.input(0)
      y = f.input(1)
      
      if max_fwd>0:
        dx = f.fwdSeed(0)
        dy = f.fwdSeed(1)
      
      z0 = 3*y
      if max_fwd>0: dz0 = 3*dy
      
      z1 = x+z0
      if max_fwd>0: dz1 = dx+dz0
      
      z2 = sin(z1)
      if max_fwd>0: dz2 = cos(z1)*dz1
      
      if max_adj>0:
        # Backwards sweep
        bx = 0
        by = 0
        bz1 = 0
        bz0 = 0
        
        bz2 = f.adjSeed(0)
        bz1 += bz2*cos(z1)
        bx+= bz1;bz0+= bz1
        by+= 3*bz0
        f.setAdjSens(bx,0)
        f.setAdjSens(by,1)
      
      f.setOutput(z2)
      if max_fwd>0: f.setFwdSens(dz2)
      
    Fun = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)] )
    Fun.init()
    with self.assertRaises(Exception):
      Fun.jacobian()
    
  def test_CustomFunction(self):
  
    x = msym("x")
    y = msym("y")
    
        
    g = MXFunction([x,y],[sin(x+3*y)])
    g.init()
    

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    def getP(max_fwd=1,max_adj=1,indirect=True):

      @pyevaluate
      def fun(f,nfwd,nadj):
        # sin(x+3*y)
        
        assert nfwd<=max_fwd
        assert nadj<=max_adj
        
        print (nfwd,nadj)
        
        x = f.input(0)
        y = f.input(1)
        
        if max_fwd>0:
          dx = f.fwdSeed(0)
          dy = f.fwdSeed(1)
        
        z0 = 3*y
        if max_fwd>0: dz0 = 3*dy
        
        z1 = x+z0
        if max_fwd>0: dz1 = dx+dz0
        
        z2 = sin(z1)
        if max_fwd>0: dz2 = cos(z1)*dz1
        
        if max_adj>0:
          # Backwards sweep
          bx = 0
          by = 0
          bz1 = 0
          bz0 = 0
          
          bz2 = f.adjSeed(0)
          bz1 += bz2*cos(z1)
          bx+= bz1;bz0+= bz1
          by+= 3*bz0
          f.setAdjSens(bx,0)
          f.setAdjSens(by,1)
        
        f.setOutput(z2)
        if max_fwd>0: f.setFwdSens(dz2)


      Fun = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)] )
      Fun.setOption("name","Fun")
      Fun.setOption("max_number_of_fwd_dir",max_fwd)
      Fun.setOption("max_number_of_adj_dir",max_adj)
      Fun.init()
      
      if not indirect: 
        Fun.setInput(0.2,0)
        Fun.setInput(0.7,1)
        return Fun

      f = MXFunction([x,y],Fun.call([x,y]))
      f.init()

      f.setInput(0.2,0)
      f.setInput(0.7,1)
      
      return f
      
    for indirect in [True,False]:
      f = getP(max_fwd=1,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
    # vector input
    
    x = msym("x",2)
    y = msym("y")
        
    g = MXFunction([x,y],[sin(x[0]+3*y)*x[1]])
    g.init()
    

    g.setInput([0.2,0.6],0)
    g.setInput(0.7,1)
    
    def getP(max_fwd=1,max_adj=1,indirect=True):

      @pyevaluate
      def fun(f,nfwd,nadj):
        # sin(x0+3*y)*x1
        
        assert nfwd<=max_fwd
        assert nadj<=max_adj
        
        print (nfwd,nadj)
        
        x0 = f.input(0)[0]
        x1 = f.input(0)[1]
        y = f.input(1)
        
        if max_fwd>0:
          dx0 = f.fwdSeed(0)[0]
          dx1 = f.fwdSeed(0)[1]
          dy = f.fwdSeed(1)
        
        z0 = 3*y
        if max_fwd>0: dz0 = 3*dy
        
        z1 = x0+z0
        if max_fwd>0: dz1 = dx0+dz0
        
        z2 = sin(z1)
        if max_fwd>0: dz2 = cos(z1)*dz1
        
        z3 = z2*x1
        if max_fwd>0: dz3 = x1*dz2 + dx1*z2
        
        if max_adj>0:
          # Backwards sweep
          bx0 = 0
          bx1 = 0
          by = 0
          
          bz2 = 0
          bz1 = 0
          bz0 = 0
          
          bz3 = f.adjSeed(0)
          bz2 += bz3*x1
          bx1 += bz3*z2
          bz1 += bz2*cos(z1)
          bx0+= bz1;bz0+= bz1
          by+= 3*bz0
          f.setAdjSens([bx0,bx1],0)
          f.setAdjSens(by,1)
        
        f.setOutput(z3)
        if max_fwd>0: f.setFwdSens(dz3)


      Fun = CustomFunction(fun, [sp_dense(2,1),sp_dense(1,1)], [sp_dense(1,1)] )
      Fun.setOption("name","Fun")
      Fun.setOption("max_number_of_fwd_dir",max_fwd)
      Fun.setOption("max_number_of_adj_dir",max_adj)
      Fun.init()

      if not indirect: 
        Fun.setInput([0.2,0.6],0)
        Fun.setInput(0.7,1)
        return Fun
        
      f = MXFunction([x,y],Fun.call([x,y]))
      f.init()

      f.setInput([0.2,0.6],0)
      f.setInput(0.7,1)
      
      return f
    
    for indirect in [True,False]:
      f = getP(max_fwd=1,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
    # vector input, vector output
    
    x = msym("x",2)
        
    g = MXFunction([x],[vertcat([x[0]**2,x[1]**2])])
    g.init()
    

    g.setInput([0.2,0.6],0)
 
    def getP(max_fwd=1,max_adj=1,indirect=True):

      @pyevaluate
      def squares(f,nfwd,nadj):
        print "Called squares with :", (nfwd,nadj)
        x = f.getInput(0)[0]
        y = f.getInput(0)[1]

        f.setOutput(f.getInput(0)**2,0)
        f.setOutput(f.getInput(0)**2,0)
        
        for i in range(nfwd):
          xdot = f.getFwdSeed(0,i)[0]
          ydot = f.getFwdSeed(0,i)[1]
          f.setFwdSens([2*x*xdot+ydot,y*xdot+x*ydot],0,i)
          
        for i in range(nadj):
          xb = f.getAdjSeed(0,i)[0]
          yb = f.getAdjSeed(0,i)[1]
          f.setAdjSens([2*x*xb+y*yb,xb+x*yb],0,i)
          
      c = CustomFunction( squares, [sp_dense(2,1)], [sp_dense(2,1)] )
      c.init()

      if not indirect: 
        c.setInput([0.2,0.6],0)
        return c
        
      f = MXFunction([x],c.call([x]))
      f.init()

      f.setInput([0.2,0.6],0)
      
      return f
    
    for indirect in [True,False]:
      f = getP(max_fwd=1,max_adj=1,indirect=indirect)
      
      with self.assertRaises(Exception):          
        self.checkfx(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
        
      with self.assertRaises(Exception): 
        self.checkfx(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      with self.assertRaises(Exception):
        self.checkfx(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
    # vector input, vector output
    
    x = msym("x",2)
        
    g = MXFunction([x],[vertcat([x[0]**2+x[1],x[0]*x[1]])])
    g.init()
    

    g.setInput([0.2,0.6],0)
 
    def getP(max_fwd=1,max_adj=1,indirect=True):

      @pyevaluate
      def squares(f,nfwd,nadj):
        print "Called squares with :", (nfwd,nadj)
        x = f.getInput(0)[0]
        y = f.getInput(0)[1]

        f.setOutput([x**2+y,x*y],0)
        
        for i in range(nfwd):
          xdot = f.getFwdSeed(0,i)[0]
          ydot = f.getFwdSeed(0,i)[1]
          f.setFwdSens([2*x*xdot+ydot,y*xdot+x*ydot],0,i)
          
        for i in range(nadj):
          xb = f.getAdjSeed(0,i)[0]
          yb = f.getAdjSeed(0,i)[1]
          f.setAdjSens([2*x*xb+y*yb,xb+x*yb],0,i)
          
      c = CustomFunction( squares, [sp_dense(2,1)], [sp_dense(2,1)] )
      c.setOption("max_number_of_fwd_dir",1)
      c.setOption("max_number_of_adj_dir",1)
      c.init()

      if not indirect: 
        c.setInput([0.2,0.6],0)
        return c
        
      f = MXFunction([x],c.call([x]))
      f.init()

      f.setInput([0.2,0.6],0)
      
      return f
    
    for indirect in [True,False]:
      f = getP(max_fwd=1,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfx(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
  def test_sparsitygenerator(self):
    x = msym("x",4)
          
      
    @sparsitygenerator
    def sp(f,iind,oind):
      return sp_dense(4,4)
      
    f = MXFunction([x],[x])
    f.init()
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.output().size(),4)
    
    f = MXFunction([x],[x])
    f.setOption("sparsity_generator",sp)
    f.init()
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.output().size(),16)
      
  def test_jacobiangenerator(self):
    x = msym("x")
    y = msym("y")
        
    g = MXFunction([x,y],[sin(x+3*y)])
    g.init()
    

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    @pyevaluate
    def fun(f,nfwd,nadj):
      # sin(x0+3*y)
      
      assert nfwd==0
      assert nadj==0
      
      x = f.input(0)
      y = f.input(1)
      
      f.setOutput(sin(x+3*y))
      
    @jacobiangenerator
    def funjac(f,iind,oind):
      # sin(x0+3*y)*x1
      print "Called jacobian with :", (iind,oind)
      
      x = ssym("x")
      y = ssym("y")
      
      if iind==0:
        f = SXFunction([x,y],[cos(x+3*y),sin(x+3*y)])
        f.init()
      elif iind==1:
        f = SXFunction([x,y],[3*cos(x+3*y),sin(x+3*y)])
        f.init()
      
      return f

    Fun = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)] )
    Fun.setOption("name","Fun")
    Fun.setOption("jacobian_generator",funjac)
    Fun.init()

    Fun.setInput(0.2,0)
    Fun.setInput(0.7,1)
    
    print Fun.input(0),Fun.input(1)
    
    print g.input(0),g.input(1)
    
    self.checkfx(Fun,g,fwd=False,adj=False,indirect=True)

    
  def test_derivative_simplifications(self):
  
    n = 1
    x = ssym("x",n)

    M = SXFunction([x],[mul((x-DMatrix(range(n))),x.T)])
    M.setOption("name","M")
    M.init()
    M.evaluate()


    P = msym("P",n,n)
    X = msym("X",n)

    MX= M.call([X])[0]

    Pf = MXFunction([X,P],[mul(MX,P)])
    Pf.setOption("name","P")
    Pf.init()

    P_P = Pf.jacobian(1)
    P_P.init()

    
    self.assertFalse("derivative" in str(P_P))
    
  def test_assert_derivatives(self):
    x = msym("x")
    
    @pyevaluate
    def dummy(f,nfwd,nadj):
      assert nfwd==0
      print f
      f.setOutput(1)

    foo = CustomFunction(dummy, [x.sparsity()], [sp_dense(1,1)] )
    foo.setOption("name","foo")
    foo.setOption("verbose",True)
    foo.setOption("max_number_of_fwd_dir",1)
    foo.init()

    y = x**2

    y = y.attachAssert(foo.call([x])[0],"P is not positive definite")

    f = MXFunction([x],[y])
    f.setOption("verbose",True)
    f.init()

    J = f.gradient()
    J.init()

    J.setInput([0.1])
    J.evaluate()
    
    print J

    self.assertFalse("derivative" in str(J))

    J = f.jacobian()
    J.init()

    J.setInput([0.1])
    J.evaluate()
    
    print J

    self.assertFalse("derivative" in str(J))
    
    H = f.hessian()
    H.init()
    
    H.setInput([0.1])
    H.evaluate()
    
if __name__ == '__main__':
    unittest.main()

