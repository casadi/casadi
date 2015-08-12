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
    fsx = SXFunction("fsx", [x,[]],[x])
    x = MX.sym("x",2)
    fmx1 = MXFunction("fmx1", [x,MX()],[x])
    fmx2 = MXFunction("fmx2", [x,[]],[x])
    
    for f in [fsx,fmx1,fmx2]:
      f.evaluate()

      X = MX.sym("X",2)
      F = f.call([X,MX()])[0]
      g = MXFunction("g", [X],[F])

      g.evaluate()
    
    x = SX.sym("x",2)
    fsx = SXFunction("fsx", [x],[x,[]])
    x = MX.sym("x",2)
    fmx1 = MXFunction("fmx1", [x],[x,MX()])
    fmx2 = MXFunction("fmx2", [x],[x,[]])
    
    for f in [fsx,fmx1,]:
      f.evaluate()

      X = MX.sym("X",2)
      F = f.call([X])[0]
      g = MXFunction("g", [X],[F])

      g.evaluate()
  
  def test_Map(self):
    self.message("Map")
    x = MX.sym("x",2)
    y = MX.sym("y")

    f = MXFunction("f", [x,y],[sin(x) + y])
        
    for mode in ["expand", "serial", "openmp"]:
      x0 = MX.sym("x0",2)
      y0 = MX.sym("y0")
      x1 = MX.sym("x1",2)
      y1 = MX.sym("y1")

      [[z0],[z1]] = f.map([[x0,y0],[x1,y1]],mode)
      
      p = MXFunction("p", [x0,y0,x1,y1],[z0,z1])
      
      n1 = DMatrix([4,5])
      N1 = 3
      n2 = DMatrix([5,7])
      N2 = 8
      

      out = p([n1,N1,n2,N2])

      self.checkarray(sin(n1)+N1,out[0],"output")
      self.checkarray(sin(n2)+N2,out[1],"output")
      
  def test_MXFunctionSeed(self):
    self.message("MXFunctionSeed")
    x1 = MX.sym("x",2)
    y1 = MX.sym("y")
    x2 = MX.sym("x",2)
    y2 = MX.sym("y")
    p= MXFunction("p", [x1,y1,x2,y2],[sin(x1) + y1,sin(x2) + y2])
    
    n1 = DMatrix([4,5])
    N1 = 3
    n2 = DMatrix([5,7])
    N2 = 8
    
    out = p([n1,N1,n2,N2])

    self.checkarray(sin(n1)+N1,out[0],"output")
    self.checkarray(sin(n2)+N2,out[1],"output")
                  
  def test_map(self):
    self.message("MX parallel call")
    x = MX.sym("x",2)
    y = MX.sym("y")

    f = MXFunction("f", [x,y],[sin(x) + y])

    #! Evaluate this function ten times in parallel
    x1 = MX.sym("x",2)
    y1 = MX.sym("y")
    x2 = MX.sym("x",2)
    y2 = MX.sym("y")
    [[F1],[F2]] = f.map([[x1,y1],[x2,y2]])
    p = MXFunction("p", [x1,y1,x2,y2],[F1,F2])
    
    n1 = DMatrix([4,5])
    N1 = 3
    n2 = DMatrix([5,7])
    N2 = 8
  
    out = p([n1,N1,n2,N2])

    self.checkarray(sin(n1)+N1,out[0],"output")
    self.checkarray(sin(n2)+N2,out[1],"output")

  def test_issue304(self):
    self.message("regression test for #304") # this code used to segfault
    x = SX.sym("x")

    f = SXFunction("f", [x],[x**2,x**3])

    X = [MX.sym("X")]

    z=f.call(X)

    g = MXFunction("g", X,[z[0]])

    g.expand([x])
  
  def test_jacobian(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)

    f = SXFunction("f", [x,y],[x**2,y,x*y[0]])

    g = f.jacobian(0,0)

    self.assertEqual(g.nIn(),f.nIn())
    self.assertEqual(g.nOut(),f.nOut()+1)

  def test_xfunction(self):
    x = SX.sym("x",3,1)
    y = SX.sym("y",2,1)
    
    f = SXFunction("f", [x,y],[x**2,y,x*y[0]])
    
    f.setInput([0.1,0.7,1.3],0)
    f.setInput([7.1,2.9],1)
    
    X = MX.sym("x",3,1)
    Y = MX.sym("y",2,1)
    
    F = MXFunction("F", [X,Y],[X**2,Y,X*Y[0]])
    
    F.setInput([0.1,0.7,1.3],0)
    F.setInput([7.1,2.9],1)
    
    self.checkfunction(f,F,sens_der=False,evals=False)
    
  
  @memory_heavy()
  def test_jacobians(self):
  
    x = SX.sym("x")
    
    self.assertEqual(jacobian(5,x).nnz(),0)
    
    
    def test(sp):
      x = SX.sym("x",sp.size2())
      sp2 = jacobian(mul(DMatrix.ones(sp),x),x).sparsity()
      self.checkarray(sp.row(),sp2.row());
      self.checkarray(sp.colind(),sp2.colind());   

    for i in range(5):
      test(Sparsity.lower(i))
      test(Sparsity.lower(i).T)
      test(Sparsity.dense(i,i))
      test(Sparsity.diag(i))
    
    for i in [63,64,65,127,128,129]:
      d = Sparsity.diag(i)
      test(d)
      
      test(d + Sparsity.rowcol([0],[5],i,i))
      
      b = Sparsity.band(i,-1) + Sparsity.band(i,1)
      test(b + Sparsity.rowcol([0],[5],i,i))
      
    m = IMatrix.ones(Sparsity.diag(129))
    m[:50,0] = 1
    m[60:,0] = 1
    m[6:9,6] = 1
    m[9,9:12] = 1
    
    sp = m[:,:120].sparsity()
    
    test(sp)
    #test(sp.T)
    
    m = IMatrix.ones(Sparsity.diag(64))
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
      test(Sparsity.lower(i))
      test(Sparsity.lower(i).T)
    
    for n in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
      for m in ([63,64,65,127,128,129] if args.run_slow else [63,64,65]):
        print (n,m)
        sp = Sparsity.dense(n,m)
        
        test(sp)
        
        random.seed(0)
        
        I = IMatrix.ones(sp)
        for i in range(n):
          for j in range(m):
            if random.random()<0.5:
              I[i,j] = 0
        I = sparsify(I)
        
        sp_holes = I.sparsity()
        
        test(sp_holes)
        
        z = IMatrix(sp_holes.size1(), sp_holes.size2())
        
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
      self.assertTrue(sp==sp.T)
      f = SXFunction("f", [x],[mul([x.T,DMatrix.ones(sp),x])])
      J = f.hessian()
      J.init()
      sp2 = J.getOutput().sparsity()
      self.checkarray(sp.row(),sp2.row())
      self.checkarray(sp.colind(),sp2.colind())
      
    A = IMatrix([[1,1,0,0,0,0],[1,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,1,0,1,1],[0,1,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()
    
    test(C)
    
    A = IMatrix([[1,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,1],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,1,0,1,0,1]])
    A = sparsify(A)
    C = A.sparsity()
      
    test(C)

    A = IMatrix([[0,0,0,0,0,0],[0,1,0,0,1,0],[0,0,1,1,0,0],[0,0,1,1,0,1],[0,1,0,0,1,0],[0,0,0,1,0,1]])
    A = sparsify(A)
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
      
      I = IMatrix.ones(sp)
      for ii in range(i):
        for jj in range(i):
          if random.random()<0.5:
            I[ii,jj] = 0
            I[jj,ii] = 0
      I = sparsify(I)
      
      sp_holes = I.sparsity()
      
      test(sp_holes)
      
      z = IMatrix(sp_holes.size1(), sp_holes.size2())
      
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
    
    f = SXFunction("f", daeIn(x=x),daeOut(ode=x))
    f.setInput([1,2])
    f.evaluate()
    a = f.getOutput()
    b = f.getOutput(0)
    c = f.getOutput("ode")
    self.checkarray(a,DMatrix([1,2]))
    self.checkarray(b,DMatrix([1,2]))
    self.checkarray(c,DMatrix([1,2]))
    f.setInput([3,4])
    f.evaluate()
    self.checkarray(a,DMatrix([1,2]))
    self.checkarray(b,DMatrix([1,2]))
    self.checkarray(c,DMatrix([1,2]))
    
  def test_customIO(self):
    
    x = SX.sym("x")
    f = SXFunction('f',[x],[x*x, x],{'output_scheme':["foo","bar"]})
    
    ret = f({"i0": 12})

    self.checkarray(DMatrix([144]),ret["foo"])
    self.checkarray(DMatrix([12]),ret["bar"])

    
    with self.assertRaises(Exception):
      f.getOutput("baz")
      
    ret = f({'i0':SX(12)})
    self.checkarray(ret["foo"],DMatrix([144]))
    self.checkarray(ret["bar"],DMatrix([12]))
    with self.assertRaises(Exception):
      self.checkarray(ret["baz"],DMatrix([12]))
     
      
  def test_unknown_options(self):
    x = SX.sym("x")
    f = SXFunction("f", [x],[x])
    
    with self.assertRaises(Exception):
      f.setOption({"fooo": False},False)
    
    f.setOption({"fooo": False},True)
    
    f.setOption({"name": "abc"},False)
    self.assertTrue(f.getOption("name")=="abc")
    f.setOption({"name": "def"},True)
    self.assertTrue(f.getOption("name")=="def")
  
  @skip("WITH_DEPRECATED_FEATURES" not in CasadiMeta.getCompilerFlags())
  def test_CustomFunctionHard(self):

    x = MX.sym("x")
    y = MX.sym("y")
    
        
    g = MXFunction("g", [x,y],[sin(x+3*y)])
    

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
          
        def getDerForward(self,f,nfwd):
          inputs = [f.getInput(i).sparsity() for i in range(f.nIn())]
          outputs = [f.getOutput(i).sparsity() for i in range(f.nOut())]
          
          sself = self

          class Der:
             def evaluate(self,xy_andseeds,z_andseeds):  sself.evaluateDerFwd(xy_andseeds,z_andseeds,nfwd)

          FunDer = PyFunction(Der(),inputs+outputs+inputs*nfwd,outputs*nfwd)
          return FunDer

        def getDerReverse(self,f,nadj):
          inputs = [f.getInput(i).sparsity() for i in range(f.nIn())]
          outputs = [f.getOutput(i).sparsity() for i in range(f.nOut())]
          
          sself = self

          class Der:
             def evaluate(self,xy_andseeds,z_andseeds):  sself.evaluateDerAdj(xy_andseeds,z_andseeds,nadj)

          FunDer = PyFunction(Der(),inputs+outputs+outputs*nadj,inputs*nadj)
          return FunDer
          
        def evaluateDerFwd(self,inputs,outputs,nfwd):
          # sin(x+3*y)
          
          num_in  =  2
          num_out =  1
          
          x = inputs[0]
          y = inputs[1]
          
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)
          
          for i in range(nfwd):
            dx = inputs[num_in + num_out + i*num_in+0]
            dy = inputs[num_in + num_out + i*num_in+1]
            
            dz0 = 3*dy
            dz1 = dx+dz0
            dz2 = cos(z1)*dz1
            
            outputs[i].set(dz2)
          
        def evaluateDerAdj(self,inputs,outputs,nadj):
          # sin(x+3*y)
          
          num_in  =  2
          num_out =  1
          
          x = inputs[0]
          y = inputs[1]
          
          z0 = 3*y
          z1 = x+z0
          z2 = sin(z1)

          
          for i in range(nadj):
            # Backwards sweep
            bx = 0
            by = 0
            bz1 = 0
            bz0 = 0
            
            bz2 = inputs[num_in + num_out + i*num_out+0]
            bz1 += bz2*cos(z1)
            bx+= bz1;bz0+= bz1
            by+= 3*bz0
            outputs[num_in*i+0].set(bx)
            outputs[num_in*i+1].set(by)
          

      Fun = PyFunction(Fun(),[Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
      if max_adj and not max_fwd:
        Fun.setOption("ad_weight", 1)
      elif max_fwd and not max_adj:
        Fun.setOption("ad_weight", 0)
      Fun.init()
      
      if not indirect: 
        Fun.setInput(0.2,0)
        Fun.setInput(0.7,1)
        return Fun

      f = MXFunction("f", [x,y],Fun.call([x,y]))

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
    
  @skip("WITH_DEPRECATED_FEATURES" not in CasadiMeta.getCompilerFlags())
  def test_CustomFunction(self):
  
    x = MX.sym("x")
    y = MX.sym("y")
    
        
    g = MXFunction("g", [x,y],[sin(x+3*y)])
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

      f = MXFunction("f", [x,y],Fun.call([x,y]))

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
      if max_adj and not max_fwd:
        Fun.setOption("ad_weight", 1)
      elif max_fwd and not max_adj:
        Fun.setOption("ad_weight", 0)
      Fun.init()
      
      if not indirect: 
        Fun.setInput(0.2,0)
        Fun.setInput(0.7,1)
        return Fun

      f = MXFunction("f", [x,y],Fun.call([x,y]))

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

      f = MXFunction("f", [x,y],Fun.call([x,y]))

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
        
    g = MXFunction("g", [x,y],[sin(x[0]+3*y)*x[1]])
    

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
      if max_adj and not max_fwd:
        Fun.setOption("ad_weight", 1)
      elif max_fwd and not max_adj:
        Fun.setOption("ad_weight", 0)
      Fun.init()

      if not indirect: 
        Fun.setInput([0.2,0.6],0)
        Fun.setInput(0.7,1)
        return Fun
        
      f = MXFunction("f", [x,y],Fun.call([x,y]))

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
        
    g = MXFunction("g", [x],[vertcat([x[0]**2+x[1],x[0]*x[1]])])
    

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
      if max_adj and not max_fwd:
        c.setOption("ad_weight", 1)
      elif max_fwd and not max_adj:
        c.setOption("ad_weight", 0)
      c.init()

      if not indirect: 
        c.setInput([0.2,0.6],0)
        return c
        
      f = MXFunction("f", [x],c.call([x]))

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
          
    f = MXFunction("f", [x],[x])
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.getOutput().nnz(),4)
    
    f = MXFunction("f", [x],[x])
    f.setJacSparsity(Sparsity.dense(4,4),0,0,True)
    
    J = f.jacobian()
    J.init()
    J.evaluate()
    
    self.assertEqual(J.getOutput().nnz(),16)
      
  def test_setjacobian(self):
    x = MX.sym("x")
    y = MX.sym("y")
        
    g = MXFunction("g", [x,y],[sin(x+3*y)])

    g.setInput(0.2,0)
    g.setInput(0.7,1)
    
    @pyevaluate
    def fun(f):
      # sin(x0+3*y)

      x = f.getInput(0)
      y = f.getInput(1)
      
      f.setOutput(sin(x+3*y))
      
    # Form Jacobians: sin(x0+3*y)*x1
    x = SX.sym("x")
    y = SX.sym("y")
    J = SXFunction("my_J", [x,y],[horzcat((cos(x+3*y),3*cos(x+3*y))),sin(x+3*y)])
    
    Fun = CustomFunction("Fun", fun, [Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)] )
    Fun.setFullJacobian(J)

    Fun.setInput(0.2,0)
    Fun.setInput(0.7,1)
    
    print Fun.getInput(0),Fun.getInput(1)
    
    print g.getInput(0),g.getInput(1)
    
    self.checkfunction(Fun,g,fwd=False,adj=False,indirect=False)

    
  def test_derivative_simplifications(self):
  
    n = 1
    x = SX.sym("x",n)

    M = SXFunction("M", [x],[mul((x-DMatrix(range(n))),x.T)])
    M.evaluate()


    P = MX.sym("P",n,n)
    X = MX.sym("X",n)

    M_X= M.call([X])[0]

    Pf = MXFunction("P", [X,P],[mul(M_X,P)])

    P_P = Pf.jacobian(1)
    P_P.init()

    
    self.assertFalse("derivative" in str(P_P))
  
  @skip("WITH_DEPRECATED_FEATURES" not in CasadiMeta.getCompilerFlags())
  def test_assert_derivatives(self):
    x = MX.sym("x")
    
    @pyevaluate
    def dummy(f):
      print f
      f.setOutput(1)

    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
        
      foo = CustomFunction(dummy, [x.sparsity()], [Sparsity.dense(1,1)] )
      foo.setOption("name","foo")
      foo.setOption("verbose",True)
      foo.init()

    # Jacobian for derivative information
    def dummy_jac(f):
      f.setOutput(1,1)

    import warnings

    with warnings.catch_warnings():
      warnings.filterwarnings("ignore",category=DeprecationWarning)
        
      foo_jac = CustomFunction(dummy_jac, [x.sparsity()], [Sparsity(1,1),Sparsity.dense(1,1)] )
      foo_jac.setOption("name","foo_jac")
      foo_jac.init()
    foo.setFullJacobian(foo_jac)

    y = x**2

    y = y.attachAssert(foo.call([x])[0],"P is not positive definite")

    f = MXFunction("f", [x], [y], {"verbose":True})

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
    
  def test_map(self):
    a = SX.sym("a",1,2)
    b = SX.sym("b")
    c = sin(a)+b
    d = cos(sumCols(a)-c)
    f = SXFunction("f",[a,b],[c,d])

    random.seed(0)

    random.random(())

    r = [[ DMatrix([1,2]).T , 3],
    [ DMatrix([2,1]).T , 1.7],
    [ DMatrix([3,4.1]).T , 2.7],
    ]

    Fref = blockcat([f(e) for e in r])


    F = MXFunction("F",[],[blockcat(f.map(r))])

    self.checkarray(F([])[0],Fref)
    
    a = SX.sym("a",1,2)
    c = sin(a)
    d = cos(sumCols(a)-c)
    f = SXFunction("f",[a],[c,d])

    random.seed(0)

    random.random(())

    r = [[ DMatrix([1,2]).T ],
    [ DMatrix([2,1]).T ],
    [ DMatrix([3,4.1]).T],
    ]

    Fref = blockcat([f(e) for e in r])


    F = MXFunction("F",[],[blockcat(f.map(r))])

    self.checkarray(F([])[0],Fref)

    a = SX.sym("a",1,2)
    b = SX.sym("b")
    c = sin(a)+b
    d = cos(sumCols(a)-c)
    f = SXFunction("f",[a,b],[c])

    random.seed(0)

    random.random(())

    r = [[ DMatrix([1,2]).T , 3],
    [ DMatrix([2,1]).T , 1.7],
    [ DMatrix([3,4.1]).T , 2.7],
    ]

    Fref = blockcat([f(e) for e in r])


    F = MXFunction("F",[],[blockcat(f.map(r))])

    self.checkarray(F([])[0],Fref)
    

  def test_simple_scheme_call(self):

    x = SX.sym("x")

    f = SXFunction("f", daeIn(x=x),[x**2])

    self.checkarray(f(x=0.3)['o0'],DMatrix(0.09))

  def test_issue1464(self):
    n = 6
    x = SX.sym("x",n)
    u = SX.sym("u")


    N = 9

    rk4 = SXFunction("f",[x,u],[x+u])
    rk4.init()

    for XX,XFunction in [(SX,SXFunction),(MX,MXFunction)]:

      g = []
      g2 = []


      V = XX.sym("V",(N+1)*n+N)
      VX,VU = vertsplit(V,[0,(N+1)*n,(N+1)*n+N])

      VXk = vertsplit(VX,n)
      VUk = vertsplit(VU,1)

      for k in range(N):
          
          [xf] = rk4([VXk[k],VUk[k]])

          xfp = vertsplit(xf,n/2)
          vp = vertsplit(VXk[k+1],n/2)

          g.append(xfp[0] - vp[0])
          g.append(xfp[1] - vp[1])

          g2.append(xf-VXk[k+1])

      for i in range(2):
        f = XFunction("nlp",[V],[vertcat(g)])
        f.setOption("ad_weight_sp",i)
        f.init()

        assert f.jacSparsity().nnz()==162

        f2 = XFunction("nlp",[V],[vertcat(g2)])
        f2.setOption("ad_weight_sp",i)
        f2.init()

        assert f2.jacSparsity().nnz()==162

  def test_callback(self):
    class mycallback(Callback2):
      def __call__(self,argin):
        return [argin[0]**2]

    c = mycallback()
    foo = c.create()
    
    x = MX.sym('x')
    y = foo([x])

    f = MXFunction("f",[x],y)
    J = f.jacobian()

    out = f([5])
    
    self.checkarray(out[0],25)

    out = J([8])
    
    self.checkarray(out[0],16,digits=6)

    class mycallback2(Callback2):
      pass

    c2 = mycallback2()
    foo = c2.create()
    
    x = MX.sym('x')
    y = foo([x])

    f = MXFunction("f",[x],y)
    J = f.jacobian()

    out = f([5])
    
    self.checkarray(out[0],10)

    out = J([8])
    
    self.checkarray(out[0],2,digits=6)

  @known_bug()
  def test_callback_errors(self):
    class mycallback(Callback2):
      def __call__(self,argin):
        raise Exception("foobar")

    c = mycallback()
    foo = c.create()
    
    x = MX.sym('x')
    y = foo([x])

    f = MXFunction("f",[x],y)

    try:
      f([3])
    except Exception as e:
      self.assertTrue("foobar" in str(e))


  def test_callback_derivatives(self):

    class mydergen(DerivativeGenerator2):
      def __init__(self,fwd=True):
        DerivativeGenerator2.__init__(self)
        self.fwd = fwd

      def __call__(self,fcn,ndir):
        # Obtain the symbols for nominal inputs/outputs
        nominal_in  = fcn.symbolicInput()
        nominal_out = fcn.symbolicOutput()

        # A growing list of inputs to the returned derivative function
        der_ins = nominal_in + nominal_out

        # A growing list of outputs to the returned derivative function
        der_outs = []

        [A] = nominal_in
        [U,S,V] = nominal_out

        constr = veccat([
          mul([U,casadi.diag(S),V])-A,
          casadi.tril(mul(U.T,U)-DMatrix.eye(mul(U.T,U).shape[0])).nz[:],
          casadi.tril(mul(V.T,V)-DMatrix.eye(mul(V.T,V).shape[0])).nz[:],
        ])

        USV = veccat(nominal_out)

        f = MXFunction("f",[USV,A],[constr])

        #impl = ImplicitFunction("impl","nlp.ipopt",f,{"nlp_solver_options" : {"print_level":0,"print_time":False}})
        impl = ImplicitFunction("impl","newton",f,{"linear_solver": "csparse"})

        if self.fwd:
          fd = impl.derivative(ndir,0)

          seeds = [ fcn.symbolicInput()[0] for i in range(ndir)]

          der_ins+=seeds

          allseeds = []
          for s in seeds:
            allseeds += [0,s]

          out = fd([USV,A] + allseeds)

          for s in out[1:]:
            [du,ds,dv]=vertsplit(s,[0,U.nnz(),U.nnz()+S.nnz(),U.nnz()+S.nnz()+V.nnz()])
            du = du.reshape(U.shape)
            ds = casadi.reshape(ds,S.shape)
            dv = casadi.reshape(dv,V.shape)

            der_outs+= [du,ds,dv]
          
        else:
          bd = impl.derivative(0,ndir)
          seeds = [ fcn.symbolicOutput() for j in range(ndir)]
          for s in seeds:
            der_ins+=s
          seedsflat = [veccat(s) for s in seeds]

          out = bd([USV,A] + seedsflat)
          

          der_outs +=out[1:][1::2]

        ret = MXFunction("my_derivative", der_ins, der_outs)
        return ret

    myd = mydergen()


    class mycallback(Callback2):
      def __init__(self,n,m, fd=True):
        Callback2.__init__(self)
        self.n = n
        self.m = m
        self.k = min(n,m)
        self.fwd = mydergen(True)
        self.adj = mydergen(False)
        self.fd = fd

      def nOut(self):
        return 3

      def inputShape(self,i):
        return (self.n,self.m)    

      def outputShape(self,i):
        if i==0:
          return (self.n, self.k) 
        elif i==1:
          return (self.k, 1)
        else:
          return (self.k, self.m)

      def __call__(self,argin):
        u,s,v = numpy.linalg.svd(argin[0],full_matrices=False)
        return [u,s,v]

      def options(self):
        return {} if self.fd else {"custom_forward": self.fwd.create(), "custom_reverse": self.adj.create()}

    n = 3
    m = 3
    c = mycallback(n,m,fd=True)
    foo = c.create()

    x = DMatrix(np.random.random((n,m)))
    X = MX.sym('x',n,m)
    Y = foo([X])

    f = MXFunction("f",[X],Y)
    Js = [f.jacobian(0,i) for i in range(3)]
    Js_ = [J([x]) for J in Js]

    u,s,v = numpy.linalg.svd(x,full_matrices=False)

    for j in Js_:
      self.checkarray(u,j[1])
      self.checkarray(s,j[2])
      self.checkarray(v,j[3])

    Js_alts = []
    for w in [0,1]:
      c = mycallback(n,m,fd=False)
      foo = c.create()

      Y = foo([X])

      f = MXFunction("f",[X],Y,{"ad_weight": w})

      J = f.jacobian(0,1)
      Js = [f.jacobian(0,i) for i in range(3)]
      Js_alt = [J([x]) for J in Js]
      Js_alts.append(Js_alt)
      for j, j_alt in zip(Js_,Js_alt):
        for i,i_alt in zip(j,j_alt):
          self.checkarray(i,i_alt,digits=5)
 
    for j, j_alt in zip(Js_alts[0],Js_alts[1]):
      for i,i_alt in zip(j,j_alt):
        self.checkarray(i,i_alt)   

  @memory_heavy()
  def test_map_node(self):
    x = SX.sym("x")
    y = SX.sym("y",2)
    z = SX.sym("z",2,2)
    v = SX.sym("z",Sparsity.upper(3))

    fun = SXFunction("f",[x,y,z,v],[mul(z,y)+x,sin(y*x).T,v/x])

    n = 2

    X = [MX.sym("x") for i in range(n)]
    Y = [MX.sym("y",2) for i in range(n)]
    Z = [MX.sym("z",2,2) for i in range(n)]
    V = [MX.sym("z",Sparsity.upper(3)) for i in range(n)]

    for parallelization in ["serial","expand","openmp"]:
      res = fun.map(zip(X,Y,Z,V),parallelization)


      flatres = []
      for r in res:
        flatres+= map(sin,r)
      F = MXFunction("F",X+Y+Z+V,flatres)

      flatresref = []
      for r in zip(X,Y,Z,V):
        flatresref+=map(sin,fun(r))

      Fref = MXFunction("F",X+Y+Z+V,flatresref)
      
      np.random.seed(0)
      X_ = [ DMatrix(i.sparsity(),np.random.random(i.nnz())) for i in X ] 
      Y_ = [ DMatrix(i.sparsity(),np.random.random(i.nnz())) for i in Y ] 
      Z_ = [ DMatrix(i.sparsity(),np.random.random(i.nnz())) for i in Z ] 
      V_ = [ DMatrix(i.sparsity(),np.random.random(i.nnz())) for i in V ] 

      for f in [F, F.expand()]:
        for i,e in enumerate(X_+Y_+Z_+V_):
          f.setInput(e,i)
          Fref.setInput(e,i)

        f.evaluate()
        Fref.evaluate()
        
        self.checkfunction(f,Fref)

  def test_issue1522(self):
    V = MX.sym("X",2)

    x =  V[0]
    y =  V[1]

    obj = (x-(x+y))**2

    nlp = MXFunction("nlp",nlpIn(x=V),nlpOut(f=obj))

    self.assertTrue(nlp.hessian(0,0).outputSparsity().issymmetric())

    V = MX.sym("X",6)

    xs =      [ V[0:2], V[2:4] ]
    travels = [ V[4],   V[5]   ]

    dist = 0

    for j in range(2):
      dist+=sumRows((xs[0]-(xs[j]+travels[j]))**2)

    nlp = MXFunction("nlp",nlpIn(x=V),nlpOut(f=-dist))

    hs = []
    for n in [nlp,SXFunction(nlp)]:
        H = n.derivative(0,1).jacobian(0,2,False,True)

        h = H(der_x=1,adj0_f=1)["jac"]
        hs.append(h)
    self.checkarray(*hs)

  def test_repmatnode(self):
    x = MX.sym("x",2)

    y = sin(repmat(x**2,1,3))

    z = MX.sym("y",2,2)

    F = MXFunction("f",[x,z],[sumCols(sumRows(y))])

    x = SX.sym("x",2)

    y = sin(repmat(x**2,1,3))
    z = SX.sym("y",2,2)

    Fref = SXFunction("f",[x,z],[sumCols(sumRows(y))])
    
    x0 = DMatrix([1,7])
    x1 = DMatrix([[3,0],[2,4]])
    F.setInput(x0)
    Fref.setInput(x0)
    F.setInput(x1,1)
    Fref.setInput(x1,1)

    self.check_codegen(F)
    self.checkfunction(F,Fref)

  def test_repsumnode(self):

    x = MX.sym("x",2)
    z = MX.sym("y",2,2)

    F = MXFunction("f",[x,z],[sin(repsum((x**2).T,1,2)),(cos(x**2)*2*x).T])

    x = SX.sym("x",2)
    z = SX.sym("y",2,2)


    Fref = SXFunction("f",[x,z],[sin(repsum((x**2).T,1,2)),(cos(x**2)*2*x).T])

    x0 = DMatrix([1,7])
    x1 = DMatrix([[3,0],[2,4]])
    F.setInput(x0)
    Fref.setInput(x0)
    F.setInput(x1,1)
    Fref.setInput(x1,1)

    self.check_codegen(F)

    self.checkfunction(F,Fref)

if __name__ == '__main__':
    unittest.main()

