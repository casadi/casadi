#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
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

class Functiontests(casadiTestCase):

  def test_call_empty(self):
    x = SX.sym("x",2)
    fsx = SXFunction([x,[]],[x])
    x = MX.sym("x",2)
    fmx1 = MXFunction([x,MX()],[x])
    fmx2 = MXFunction([x,[]],[x])
    
    for f in [fsx,fmx1,fmx2]:
      f.init()
      f.evaluate()

      X = MX.sym("X",2)
      F = f.call([X,X])[0]
      g = MXFunction([X],[F])
      g.init()

      g.evaluate()
    
    x = SX.sym("x",2)
    fsx = SXFunction([x],[x,[]])
    x = MX.sym("x",2)
    fmx1 = MXFunction([x],[x,MX()])
    fmx2 = MXFunction([x],[x,[]])
    
    for f in [fsx,fmx1,]:
      f.init()
      f.evaluate()

      X = MX.sym("X",2)
      F = f.call([X])[0]
      g = MXFunction([X],[F])
      g.init()

      g.evaluate()
  
  def test_Parallelizer(self):
    self.message("Parallelizer")
    x = MX.sym("x",2)
    y = MX.sym("y")

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
      
      p.evaluate()

      self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
      self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
      
  def test_MXFunctionSeed(self):
    self.message("MXFunctionSeed")
    x1 = MX.sym("x",2)
    y1 = MX.sym("y")
    x2 = MX.sym("x",2)
    y2 = MX.sym("y")
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
        
    p.evaluate()

    self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
    self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
        
  def test_Parallelizer2(self):
    self.message("Parallelizer")
    x = MX.sym("x",2)
    y = MX.sym("y")

    f = MXFunction([x,y],[sin(x) + y])
    f.init()
    
    #! Evaluate this function ten times in parallel
    pp = Parallelizer([f]*2)
    for mode in ["serial","openmp"]:
      pp.setOption("parallelization",mode)
      pp.init()
      
      x1 = MX.sym("x",2)
      y1 = MX.sym("y")
      x2 = MX.sym("x",2)
      y2 = MX.sym("y")
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
      
      p.evaluate()

      self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
      self.checkarray(sin(n2)+N2,p.getOutput(1),"output")
          
  def test_ParallelizerMXCall(self):
    self.message("MX parallel call")
    x = MX.sym("x",2)
    y = MX.sym("y")

    f = MXFunction([x,y],[sin(x) + y])
    f.init()

    #! Evaluate this function ten times in parallel
    x1 = MX.sym("x",2)
    y1 = MX.sym("y")
    x2 = MX.sym("x",2)
    y2 = MX.sym("y")
    [[F1],[F2]] = f.callParallel([[x1,y1],[x2,y2]])
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
    
    p.evaluate()

    self.checkarray(sin(n1)+N1,p.getOutput(0),"output")
    self.checkarray(sin(n2)+N2,p.getOutput(1),"output")

  def test_set_wrong(self):
    self.message("setter, wrong sparsity")
    x = SXElement.sym("x")

    f = SXFunction([x],[x])
    f.init()
    A = DMatrix.sparse(1,1)
    #self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0)) # This is now o.k. syntax
    B = DMatrix.repmat(2,1,2)
    #self.assertRaises(RuntimeError,lambda : f.getFwdSeed(A,0)) # This is now o.k. syntax
    
  def test_issue304(self):
    self.message("regression test for #304") # this code used to segfault
    x = SXElement.sym("x")

    f = SXFunction([x],[x**2,x**3])
    f.init()

    X = [MX.sym("X")]

    z=f.call(X)

    g = MXFunction(X,[z[0]])
    g.init()

    g.expand([x])
  
  def test_jacobian(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)

    f = SXFunction([x,y],[x**2,y,x*y[0]])
    f.init()

    g = f.jacobian(0,0)

    self.assertEqual(g.getNumInputs(),f.getNumInputs())
    self.assertEqual(g.getNumOutputs(),f.getNumOutputs()+1)

  def test_xfunction(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)
    
    f = SXFunction([x,y],[x**2,y,x*y[0]])
    f.init()
    
    f.setInput([0.1,0.7,1.3],0)
    f.setInput([7.1,2.9],1)
    
    X = MX.sym("x",3,1)
    Y = MX.sym("y",2,1)
    
    F = MXFunction([X,Y],[X**2,Y,X*Y[0]])
    F.init()
    
    F.setInput([0.1,0.7,1.3],0)
    F.setInput([7.1,2.9],1)
    
    self.checkfunction(f,F,sens_der=False,evals=False)
    
  
  @memory_heavy()
  def test_jacobians(self):
  
    x = SX.sym("x")
    
    self.assertEqual(jacobian(5,x).size(),0)
    
    
    def test(sp):
      x = SX.sym("x",sp.size2())
      sp2 = jacobian(mul(DMatrix(sp,1),x),x).sparsity()
      self.checkarray(sp.row(),sp2.row());
      self.checkarray(sp.colind(),sp2.colind());   

    for i in range(5):
      test(Sparsity.tril(i))
      test(Sparsity.tril(i).T)
      test(Sparsity.dense(i,i))
      test(Sparsity.diag(i))
    
    for i in [63,64,65,127,128,129]:
      d = Sparsity.diag(i)
      test(d)
      
      test(d + Sparsity.rowcol([0],[5],i,i))
      
      b = Sparsity.band(i,-1) + Sparsity.band(i,1)
      test(b + Sparsity.rowcol([0],[5],i,i))
      
    m = IMatrix(Sparsity.diag(129),1)
    m[:50,0] = 1
    m[60:,0] = 1
    m[6:9,6] = 1
    m[9,9:12] = 1
    
    sp = m[:,:120].sparsity()
    
    test(sp)
    #test(sp.T)
    
    m = IMatrix(Sparsity.diag(64),1)
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
      test(Sparsity.tril(i))
      test(Sparsity.tril(i).T)
    
    for n in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
      for m in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
        print (n,m)
        sp = Sparsity.dense(n,m)
        
        test(sp)
        
        random.seed(0)
        
        I = IMatrix(sp,1)
        for i in range(n):
          for j in range(m):
            if random.random()<0.5:
              I[i,j] = 0
        I = sparse(I)
        
        sp_holes = I.sparsity()
        
        test(sp_holes)
        
        z = IMatrix.sparse(sp_holes.shape)
        
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
      x = SX.sym("x",sp.size2())
      self.assertTrue(sp==sp.transpose())
      f = SXFunction([x],[mul([x.T,DMatrix(sp,1),x])])
      f.init()
      J = f.hessian()
      J.init()
      sp2 = J.output().sparsity()
      self.checkarray(sp.row(),sp2.row())
      self.checkarray(sp.colind(),sp2.colind())
      
    A = IMatrix([[1,1,0,0,0,0],[1,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparse(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparse(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparse(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparse(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,0,0,1,0,1]])
    A = sparse(A)
    C = A.sparsity()
      
    test(C)
    
    
    for i in [63,64,65,100,127,128,129]:
      d = Sparsity.diag(i)
      test(d)
      
      test(d + Sparsity.rowcol([0],[5],i,i) + Sparsity.rowcol([5],[0],i,i))
      
      b = Sparsity.band(i,-1) + Sparsity.band(i,1)
      test(b)
      test(b + Sparsity.rowcol([0],[5],i,i) + Sparsity.rowcol([5],[0],i,i))
      
      d = Sparsity.dense(i,i)
      test(d)
      
      d = Sparsity.diag(i) + Sparsity.triplet(i,i,[0]*i,range(i))+Sparsity.triplet(i,i,range(i),[0]*i)
      test(d)


      sp = Sparsity.dense(i,i)
        
      random.seed(0)
      
      I = IMatrix(sp,1)
      for ii in range(i):
        for jj in range(i):
          if random.random()<0.5:
            I[ii,jj] = 0
            I[jj,ii] = 0
      I = sparse(I)
      
      sp_holes = I.sparsity()
      
      test(sp_holes)
      
      z = IMatrix.sparse(sp_holes.shape)
      
      R = 5
      v = []
      for r in range(R):
        h = [z]*5
        h[r] = I
        v.append(horzcat(h))
      d = vertcat(v)
      
      test(d.sparsity())
      
  def test_getOutput(self):
    x = SX.sym("x",2)
    
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

    x = SX.sym("x")
    
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
      
    ret = f.call([SX(12)])
    self.checkarray(myOut(ret,"foo")[0],DMatrix([144]))
    self.checkarray(myOut(ret,"bar")[0],DMatrix([12]))
    with self.assertRaises(Exception):
      self.checkarray(myOut(ret,"baz")[0],DMatrix([12]))
     
      
  def test_unknown_options(self):
    x = SX.sym("x")
    f = SXFunction([x],[x])
    f.init()
    
    with self.assertRaises(Exception):
      f.setOption({"fooo": False},False)
    
    f.setOption({"fooo": False},True)
    
    f.setOption({"name": "abc"},False)
    self.assertTrue(f.getOption("name")=="abc")
    f.setOption({"name": "def"},True)
    self.assertTrue(f.getOption("name")=="def")
    
  def test_CustomFunctionHard(self):
  
    x = MX.sym("x")
    y = MX.sym("y")
    
        
    g = MXFunction([x,y],[sin(x+3*y)])
    g.init()
    

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    def getP(max_fwd=1,max_adj=1,indirect=True):

      class Fun:
        # sin(x+3*y)
        
        def evaluate(self,(x,y),(z,)):
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          z.set(z2)
          
        def getDerivative(self,f,nfwd,nadj):
          inputs = [f.input(i).sparsity() for i in range(f.getNumInputs())]
          outputs = [f.output(i).sparsity() for i in range(f.getNumOutputs())]
          
          sself = self

          class Der:
             def evaluate(self,xy_andseeds,z_andseeds):  sself.evaluateDer(xy_andseeds,z_andseeds,nfwd,nadj)

          FunDer = PyFunction(Der(),inputs+inputs*nfwd+outputs*nadj,outputs*(nfwd+1)+inputs*nadj)
          return FunDer
          
        def evaluateDer(self,inputs,outputs,nfwd,nadj):
          # sin(x+3*y)
          
          num_in  =  2
          num_out =  1
          
          x = inputs[0]
          y = inputs[1]
          
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          outputs[0].set(z2)
          
          for i in range(nfwd):
            dx = inputs[num_in + i*num_in+0]
            dy = inputs[num_in + i*num_in+1]
            
            dz0 = 3*dy
            dz1 = dx+dz0
            dz2 = cos(z1)*dz1
            
            outputs[num_out + i].set(dz2)
          
          for i in range(nadj):
            # Backwards sweep
            bx = 0
            by = 0
            bz1 = 0
            bz0 = 0
            
            bz2 = inputs[num_in + nfwd*num_in+i*num_out+0]
            bz1 += bz2*cos(z1)
            bx+= bz1;bz0+= bz1
            by+= 3*bz0
            outputs[num_out + nfwd*num_out + num_in*i+0].set(bx)
            outputs[num_out + nfwd*num_out + num_in*i+1].set(by)
          

      Fun = PyFunction(Fun(),[Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
      if max_fwd and max_adj:
        Fun.setOption("ad_mode","automatic")
      elif max_adj:
        Fun.setOption("ad_mode","reverse")
      elif max_fwd:
        Fun.setOption("ad_mode","forward")
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
                
      self.checkfunction(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
    
  def test_CustomFunction(self):
  
    x = MX.sym("x")
    y = MX.sym("y")
    
        
    g = MXFunction([x,y],[sin(x+3*y)])
    g.init()
    

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    # Simple syntax
    def getP(indirect=True):
    
      @pyfunction([Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
      def Fun((x,y)):
        # sin(x+3*y)
        
        z0 = 3*y
        z1 = x+z0
        z2 = sin(z1)
        return [z2]
          
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
      f = getP(indirect=indirect)
      self.checkfunction(f,g,sens_der=False,jacobian=False,gradient=False,hessian=False,evals=False)
    
      with self.assertRaises(Exception):          
        f.gradient()

      with self.assertRaises(Exception): 
        f.jacobian()
    
      with self.assertRaises(Exception):
        f.derivative()
        
    def getP(max_fwd=1,max_adj=1,indirect=True):

      class Fun:
        # sin(x+3*y)
        
        def evaluate(self,(x,y),(z,)):
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          z.set(z2)

        def fwd(self,(x,y),(z,),seeds,sens):
          assert(max_fwd)
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          z.set(z2)
          
          for ((dx,dy),(dz,)) in zip(seeds,sens):
            dz0 = 3*dy
            dz1 = dx+dz0
            dz2 = cos(z1)*dz1
            dz.set(dz2)
        
        def adj(self,(x,y),(z,),seeds,sens):
          assert(max_adj)
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          z.set(z2)
          
          for ((z_bar,),(x_bar,y_bar)) in zip(seeds,sens):
            bx = 0
            by = 0
            bz1 = 0
            bz0 = 0
            
            bz2 = z_bar
            bz1 += bz2*cos(z1)
            bx+= bz1;bz0+= bz1
            by+= 3*bz0
            x_bar.set(bx)
            y_bar.set(by)

      Fun = PyFunction(Fun(),[Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
      if max_fwd and max_adj:
        Fun.setOption("ad_mode","automatic")
      elif max_adj:
        Fun.setOption("ad_mode","reverse")
      elif max_fwd:
        Fun.setOption("ad_mode","forward")
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
                
      self.checkfunction(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,fwd=False,evals=1)

    def getP(max_fwd=1,max_adj=1,indirect=True):

      class Fun:
        # sin(x+3*y)
        
        def evaluate(self,(x,y),(z,)):
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          z.set(z2)
          
        if max_fwd:
          def fwd(self,(x,y),(z,),seeds,sens):
            z0 = 3*y
            z1 = x+z0
            z2 = sin(z1)
            z.set(z2)
            
            for ((dx,dy),(dz,)) in zip(seeds,sens):
              dz0 = 3*dy
              dz1 = dx+dz0
              dz2 = cos(z1)*dz1
              dz.set(dz2)
        
        if max_adj:
          def adj(self,(x,y),(z,),seeds,sens):
            z0 = 3*y
            z1 = x+z0
            z2 = sin(z1)
            z.set(z2)
            
            for ((z_bar,),(x_bar,y_bar)) in zip(seeds,sens):
              bx = 0
              by = 0
              bz1 = 0
              bz0 = 0
              
              bz2 = z_bar
              bz1 += bz2*cos(z1)
              bx+= bz1;bz0+= bz1
              by+= 3*bz0
              x_bar.set(bx)
              y_bar.set(by)

      Fun = PyFunction(Fun(),[Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
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
                
      self.checkfunction(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
    # vector input
    
    x = MX.sym("x",2)
    y = MX.sym("y")
        
    g = MXFunction([x,y],[sin(x[0]+3*y)*x[1]])
    g.init()
    

    g.setInput([0.2,0.6],0)
    g.setInput(0.7,1)
    
    def getP(max_fwd=1,max_adj=1,indirect=True):

      class Fun:
        
        def evaluate(self,(x,y),(z,)):
          # sin(x0+3*y)*x1

          x0 = x[0]
          x1 = x[1]
          
          z0 = 3*y
          z1 = x0+z0
          z2 = sin(z1)
          z3 = z2*x1
          
          z.set(z3)
        
        def fwd(self,(x,y),(z,),seeds,sens):
          assert(max_fwd)
          x0 = x[0]
          x1 = x[1]
          
          z0 = 3*y
          z1 = x0+z0
          z2 = sin(z1)
          z3 = z2*x1
          
          z.set(z3)
          
          for ((dx,dy),(dz,)) in zip(seeds,sens):
            dx0=dx[0]
            dx1=dx[1]

            dz0 = 3*dy
            dz1 = dx0+dz0
            dz2 = cos(z1)*dz1
            dz3 = x1*dz2 + dx1*z2
          
            dz.set(dz3)
        def adj(self,(x,y),(z,),seeds,sens):
          assert(max_adj)
          x0 = x[0]
          x1 = x[1]
          
          z0 = 3*y
          z1 = x0+z0
          z2 = sin(z1)
          z3 = z2*x1
          
          z.set(z3)
          
          for ((z_bar,),(x_bar,y_bar)) in zip(seeds,sens):
            # Backwards sweep
            bx0 = 0
            bx1 = 0
            by = 0
            
            bz2 = 0
            bz1 = 0
            bz0 = 0
            
            bz3 = z_bar
            bz2 += bz3*x1
            bx1 += bz3*z2
            bz1 += bz2*cos(z1)
            bx0+= bz1;bz0+= bz1
            by+= 3*bz0
            x_bar.set([bx0,bx1])
            y_bar.set(by)

      Fun = PyFunction(Fun(),[Sparsity.dense(2,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
      if max_fwd and max_adj:
        Fun.setOption("ad_mode","automatic")
      elif max_adj:
        Fun.setOption("ad_mode","reverse")
      elif max_fwd:
        Fun.setOption("ad_mode","forward")
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
                
      self.checkfunction(f,g,sens_der=False,hessian=False,evals=1)

      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
      
    # vector input, vector output
    
    x = MX.sym("x",2)
        
    g = MXFunction([x],[vertcat([x[0]**2+x[1],x[0]*x[1]])])
    g.init()
    

    g.setInput([0.2,0.6],0)
 
    def getP(max_fwd=1,max_adj=1,indirect=True):

      
      class Squares:

         def evaluate(self,(X,),(Y,)):
            x = X[0]
            y = X[1]
            Y.set([x**2+y,x*y])
          
         def fwd(self,(X,),(Y,),seeds,sens):
            assert(max_fwd)
            x = X[0]
            y = X[1]
            Y.set([x**2+y,x*y])
            for ((Xdot,),(Zdot,)) in zip(seeds,sens):
              xdot = Xdot[0]
              ydot = Xdot[1]
              Zdot.set([2*x*xdot+ydot,y*xdot+x*ydot])
            
         def adj(self,(X,),(Y,),seeds,sens):
            assert(max_adj)
            x = X[0]
            y = X[1]
            Y.set([x**2+y,x*y])
            for ((Y_bar,),(X_bar,)) in zip(seeds,sens):
              xb = Y_bar[0]
              yb = Y_bar[1]
              X_bar.set([2*x*xb+y*yb,xb+x*yb])
          
      c = PyFunction(Squares(),[Sparsity.dense(2,1)], [Sparsity.dense(2,1)])
      if max_fwd and max_adj:
        c.setOption("ad_mode","automatic")
      elif max_adj:
        c.setOption("ad_mode","reverse")
      elif max_fwd:
        c.setOption("ad_mode","forward")
      c.init()

      if not indirect: 
        c.setInput([0.2,0.6],0)
        return c
        
      f = MXFunction([x],c.call([x]))
      f.init()

      f.setInput([0.2,0.6],0)
      
      return f
    
    for indirect in [False]:
      f = getP(max_fwd=1,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,evals=1)

      print f
      f = getP(max_fwd=1,max_adj=0,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,adj=False,evals=1)

      f = getP(max_fwd=0,max_adj=1,indirect=indirect)
                
      self.checkfunction(f,g,sens_der=False,hessian=False,fwd=False,evals=1)
         
  def test_setjacsparsity(self):
    x = MX.sym("x",4)
          
    f = MXFunction([x],[x])
    f.init()
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.output().size(),4)
    
    f = MXFunction([x],[x])
    f.init()
    f.setJacSparsity(Sparsity.dense(4,4),0,0,True)
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.output().size(),16)
      
  def test_setjacobian(self):
    x = MX.sym("x")
    y = MX.sym("y")
        
    g = MXFunction([x,y],[sin(x+3*y)])
    g.init()
    

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    @pyevaluate
    def fun(f):
      # sin(x0+3*y)

      x = f.input(0)
      y = f.input(1)
      
      f.setOutput(sin(x+3*y))
      
    # Form Jacobians: sin(x0+3*y)*x1
    x = SX.sym("x")
    y = SX.sym("y")
    J = SXFunction([x,y],[horzcat((cos(x+3*y),3*cos(x+3*y))),sin(x+3*y)])
    J.setOption("name","my_J")
    J.init()
    
    Fun = CustomFunction(fun, [Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)] )
    Fun.setOption("name","Fun")
    Fun.init()
    Fun.setFullJacobian(J)

    Fun.setInput(0.2,0)
    Fun.setInput(0.7,1)
    
    print Fun.input(0),Fun.input(1)
    
    print g.input(0),g.input(1)
    
    self.checkfunction(Fun,g,fwd=False,adj=False,indirect=False)

    
  def test_derivative_simplifications(self):
  
    n = 1
    x = SX.sym("x",n)

    M = SXFunction([x],[mul((x-DMatrix(range(n))),x.T)])
    M.setOption("name","M")
    M.init()
    M.evaluate()


    P = MX.sym("P",n,n)
    X = MX.sym("X",n)

    M_X= M.call([X])[0]

    Pf = MXFunction([X,P],[mul(M_X,P)])
    Pf.setOption("name","P")
    Pf.init()

    P_P = Pf.jacobian(1)
    P_P.init()

    
    self.assertFalse("derivative" in str(P_P))
    
  def test_assert_derivatives(self):
    x = MX.sym("x")
    
    @pyevaluate
    def dummy(f):
      print f
      f.setOutput(1)

    foo = CustomFunction(dummy, [x.sparsity()], [Sparsity.dense(1,1)] )
    foo.setOption("name","foo")
    foo.setOption("verbose",True)
    foo.init()

    # Jacobian for derivative information
    def dummy_jac(f):
      f.setOutput(1,1)

    foo_jac = CustomFunction(dummy_jac, [x.sparsity()], [Sparsity.sparse(1,1),Sparsity.dense(1,1)] )
    foo_jac.setOption("name","foo_jac")
    foo_jac.init()
    foo.setFullJacobian(foo_jac)

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
    
  def test_simple_scheme_call(self):

    x = SX.sym("x")

    f = SXFunction(daeIn(x=x),[x**2])
    f.init()

    self.checkarray(f(x=0.3)[0],DMatrix(0.09))
    
if __name__ == '__main__':
    unittest.main()

