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
from casadi.tools import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

sdpsolvers = []
try:
  sdpsolvers.append((DSDPSolver,{}))
except:
  pass
  
class SDPtests(casadiTestCase):

  def test_memleak1(self):
    self.message("memleak1")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    
    A = DMatrix(0,3)
    
    c = DMatrix([48,-8,20])

    F = -vertcat([DMatrix([[10,4],[4,0]]),DMatrix([[0,0],[0,-8]]),DMatrix([[0,-8],[-8,-2]])])

    makeSparse(F)

    F.printMatrix()

    G = -DMatrix([[-11,0],[0,23]])

    makeSparse(G)
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))

  def test_memleak2(self):
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    
    A = DMatrix(0,3)
    
    c = DMatrix([48,-8,20])

    F = -vertcat([DMatrix([[10,4],[4,0]]),DMatrix([[0,0],[0,-8]]),DMatrix([[0,-8],[-8,-2]])])

    makeSparse(F)

    F.printMatrix()

    G = -DMatrix([[-11,0],[0,23]])

    makeSparse(G)

    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()

  def test_nodsdp(self):
    self.message("scalar")
    

    n1 = 3.1
    c = DMatrix(n1)
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=sp_dense(0,1),g=sp_dense(0,0),f=sp_dense(0,0)))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(c,"c")
      sdp.setInput(-1,"lbx")
      sdp.setInput(1,"ubx")
      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),-n1,digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),-n1,digits=5)
      self.checkarray(sdp.getOutput("x"),-1,digits=5)

  def test_simple_sdp_A(self):
    self.message("scalar")
    
    A = DMatrix(0,1)
     
    #
    # min  x
    #  x
    #      |-2   x| <=0
    #      |x   -2|
    #
    #
    #   shur: -2 < 0
    #         -2 - x^2/(-2) < 0
    #         4 - x^2 > 0
    #         x^2 < 4
    #     ->    -2 <= x <= 2 
    #  
    c = DMatrix([1])
    Fi = [DMatrix([[0,1],[1,0]])]
    F = vertcat(Fi)
    G = DMatrix([[2,0],[0,2]])
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")
      
      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("x"),DMatrix([-2]),digits=5)
      self.checkarray(sdp.getOutput("cost"),DMatrix([-2]),digits=5)


  def test_simple_sdp(self):
    self.message("scalar")
    
    A = DMatrix(0,2)
     
    #
    # min  2*x+y
    #  x,y
    #      |-x   2| <=0
    #      |2   -y|
    #   x,y>=0
    #
    #
    #   shur: -x < 0
    #         -y - 4/(-x) < 0
    #         xy - 4 > 0
    #         xy > 4
    #  
    #     <=>    border: xy - 4 == 0
    # Max(Eigenvalues({-x,2},{2,-y}))    
    #     <=>    xy >= 4
    #
    #  [-1 0;0 0] x + [0 0;0 -1] y - [0 -2; -2 0]
    c = DMatrix([2,1])
    Fi = [DMatrix([[-1,0],[0,0]]),DMatrix([[0,0],[0,-1]])]
    F = vertcat(Fi)
    G = DMatrix([[0,-2],[-2,0]])
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")
      sdp.setInput(0,"lbx")
      sdp.setInput(10,"ubx")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("x"),DMatrix([sqrt(2),2*sqrt(2)]),digits=5)
      self.checkarray(sdp.getOutput("cost"),DMatrix([2*sqrt(2)+2*sqrt(2)]),digits=5)

  def test_scalar(self):
    self.message("scalar")
    
    A = DMatrix(0,1)
     
    #
    # min  n1*x
    #  x
    #       n3*x-n2>=0
    #
    #  -> x = n2/n3
    #    
    #  1 active constraint, cost:  d(n1*x)/d(n*x-n2) = 1/[d(n3*x-n2)/d(n1*x)] = n1/n3
    n1 = 3.1
    n2 = 2.3
    n3 = 4.7
    c = DMatrix(n1)
    Fi = [DMatrix(n3)]
    F = -vertcat(Fi)
    makeSparse(F)
    G = -DMatrix(n2)
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(n1*n2/n3),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(n1*n2/n3),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix(n2/n3),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix(0),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),DMatrix(n1/n3),digits=5)

  def test_linear_equality(self):
  
    A = DMatrix(0,1)
    self.message("linear equality")
    
    #  min   n1*x
    #   x
    #
    #    n3*x-n2    >= 0   |__   n3*x == n2
    #    -(n3*x-n2) >= 0   |
    #
    # solution: x=n2/n3
    
    n3 = 1.7
    n1 = 2.1
    n2 = 1.3
    c = DMatrix([n1])
    Fi = [ blkdiag([n3,-n3])]
    G = -blkdiag([n2,-n2])
    F = -vertcat(Fi)
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(n1*n2/n3),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(n1*n2/n3),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix(n2/n3),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix.zeros(2,2),digits=5)
      
      self.checkarray(sdp.getOutput("dual")[0,0]-sdp.getOutput("dual")[1,1],DMatrix(n1/n3),digits=5)

  def test_linear_interpolation1(self):
    self.message("linear interpolation1")

    #  min    2*x0 + x1*3
    #   x0,x1
    #          x0+x1 - 1 >=0  -->  x0+x1>=1
    #            x0  >=0
    #            x1  >=0
    #                 
    #  solution: x0=1, x1=0
    
    A = DMatrix(0,2)
    
    c = DMatrix([2,3])
    Fi = [ blkdiag([1,1,0]), blkdiag([1,0,1])]
    G = -blkdiag([1,0,0])
    F = -vertcat(Fi)
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([1,0]),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix([[0,0,0],[0,1,0],[0,0,0]]),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),DMatrix([[2,0,0],[0,0,0],[0,0,1]]),digits=5)

  def test_linear_interpolation1_bounds(self):
    self.message("linear interpolation1")

    #  min    2*x0 + x1*3
    #   x0,x1
    #          x0+x1 - 1 >=0  -->  x0+x1>=1
    #            x0  >=0
    #            x1  >=0
    #                 
    #  solution: x0=1, x1=0
    
    A = DMatrix(0,2)
    
    c = DMatrix([2,3])
    Fi = [ blkdiag([1]), blkdiag([1])]
    G = -blkdiag([1])
    F = -vertcat(Fi)
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")
      sdp.setInput([0,0],"lbx")
      sdp.setInput([inf,inf],"ubx")
      
      sdp.evaluate()
      self.checkarray(sdp.getOutput("x"),DMatrix([1,0]),digits=5)
      self.checkarray(sdp.getOutput("cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix([0]),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),DMatrix([2]),digits=5)
      
      self.checkarray(sdp.getOutput("lam_x"),DMatrix([0,-1]),digits=5)

  def test_linear_interpolation1_A_boundsviol(self):

    A = DMatrix([[1,1]])
    c = DMatrix([2,3])
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=sp_dense(0,0),f=sp_dense(0,0)))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(c,"c")
      sdp.setInput(A,"a")
      sdp.setInput([0,0],"lbx")
      sdp.setInput([inf,-3],"ubx")
      sdp.setInput([1],"lba")
      sdp.setInput([inf],"uba")
      
      with self.assertRaises(Exception):
        sdp.evaluate()
        
      sdp.setInput([0,0],"lbx")
      sdp.setInput([inf,inf],"ubx")
      sdp.setInput([9],"lba")
      sdp.setInput([6],"uba")
      
      with self.assertRaises(Exception):
        sdp.evaluate()
    
  def test_linear_interpolation1_A(self):
    self.message("linear interpolation1")

    #  min    2*x0 + x1*3
    #   x0,x1
    #          x0+x1 - 1 >=0  -->  x0+x1>=1
    #            x0  >=0
    #            x1  >=0
    #                 
    #  solution: x0=1, x1=0
    
    A = DMatrix([[1,1]])
    c = DMatrix([2,3])
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=sp_dense(0,0),f=sp_dense(0,0)))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(c,"c")
      sdp.setInput(A,"a")
      sdp.setInput([0,0],"lbx")
      sdp.setInput([inf,inf],"ubx")
      sdp.setInput([1],"lba")
      sdp.setInput([inf],"uba")
      
      sdp.evaluate()
      self.checkarray(sdp.getOutput("x"),DMatrix([1,0]),digits=5)
      self.checkarray(sdp.getOutput("cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("lam_x"),DMatrix([0,-1]),digits=5)
      self.checkarray(sdp.getOutput("lam_a"),DMatrix([-2]),digits=5)
    
  def test_linear_interpolation2(self):
    self.message("linear interpolation2")

    A = DMatrix(0,2)
     
    #  min     2*x0 + 3*x1
    #   x0,x1
    #           -(x0 + x1 -1) >=0  -->  x0 + x1 <= 1
    #                 x0 >=0
    #                 x1 >=0
    #
    # solution:  x0=0 , x1=0
    c = DMatrix([2,3])
    Fi = [ blkdiag([-1,1,0]), blkdiag([-1,0,1])]
    G = -blkdiag([-1,0,0])
    F = -vertcat(Fi)
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(0),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(0),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([0,0]),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix([[1,0,0],[0,0,0],[0,0,0]]),digits=5)
      self.checkarray(sdp.getOutput("dual"),DMatrix([[0,0,0],[0,2,0],[0,0,3]]),digits=5)

  def test_linear_interpolation(self):
    self.message("linear interpolation")
    
    A = DMatrix(0,2)
    
    #  min  2*a + (1-a)*4
    #   a
    #          0  <= a <=  1
    #                 
    
    
    # Translates to: 
    #    min     2*x0 + 4*x1
    #   x0,x1
    #           x0 + x1 -1  >= 0  |__   x0 + x1 == 1
    #         -(x0 + x1 -1) >= 0  |
    #                x0     >= 0    
    #                x1     >= 0
    
    c = DMatrix([2,4])
    Fi = [ blkdiag([1,-1,1,0]), blkdiag([1,-1,0,1])]
    e = 1e-6
    G = -blkdiag([1,-(1+e),0,0])
    F = -vertcat(Fi)
    
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()
      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(2),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([1,0]),digits=5)
      self.checkarray(sdp.getOutput("p"),diag([0,0,1,0]),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),diag([2,0,0,2]),digits=2)

  def test_example1(self):
    self.message("Example1")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    c = DMatrix([48,-8,20])
    
    A = DMatrix(0,3)
    
    Fi = [-DMatrix([[10,4],[4,0]]),-DMatrix([[0,0],[0,-8]]),-DMatrix([[0,-8],[-8,-2]])]

    F = vertcat(Fi)

    makeSparse(F)

    F.printMatrix()

    G = -DMatrix([[-11,0],[0,23]])

    makeSparse(G)

    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()

      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      
      self.checkarray(sdp.getOutput("cost"),DMatrix(-41.9),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(-41.9),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([-1.1,-2.7375,-0.55]),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),DMatrix([[5.9,-1.375],[-1.375,1]]),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix.zeros(2,2),digits=5)
      
      try:
        IpoptSolver
      except:
        return
        
      V = struct_ssym([
            entry("L",shape=G.shape),
            entry("x",shape=c.size())
          ])
      L = V["L"]
      x = V["x"] 

      P = mul(L,L.T)


      g = []
      g.append(sum([-Fi[i]*x[i] for i in range(3)]) + G - P)

      nlp = SXFunction(nlpIn(x=V),nlpOut(f=mul(c.T,x),g=veccat(g)))

      sol = IpoptSolver(nlp)
      sol.init()
      sol.setInput(0,"lbg")
      sol.setInput(0,"ubg")
      sol.setInput(1,"x0")

      sol.evaluate()

      sol_ = V(sol.getOutput())
      
      self.checkarray(sol_["x"],DMatrix([-1.1,-2.7375,-0.55]),digits=5)
      

  def test_example2(self):
    self.message("Example2")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    c = DMatrix([1.1, -10, 6.6 , 19 , 4.1])

    A = DMatrix(0,5)

    G = -blkdiag([DMatrix([[-1.4,-3.2],[-3.2,-28]]),DMatrix([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0]);
    
    sp = G.sparsity()
    
    flatdata = [[0.5,5.2,5.2,-5.3,7.8,-2.4,6.0,-2.4,4.2,6.5,6.0,6.5,2.1,-4.5,-3.5],
    [1.7,7.0,7.0,-9.3,-1.9,-0.9,-1.3,-0.9,-0.8,-2.1,-1.3,-2.1,4.0,-0.2,-3.7],
 [6.3,-7.5,-7.5,-3.3,0.2,8.8,5.4,8.8,3.4,-0.4,5.4,-0.4,7.5,-3.3,-4.0],
  [-2.4,-2.5,-2.5,-2.9,3.4,-3.2,-4.5,-3.2,3.0,-4.8,-4.5,-4.8,3.6,4.8,9.7],
  [-6.5,-5.4,-5.4,-6.6,6.7,-7.2,-3.6,-7.2,7.3,-3.0,-3.6,-3.0,-1.4,6.1,-1.5]]

    F = -vertcat([DMatrix(sp,data) for data in flatdata])
    makeSparse(F)


    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()

      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      DMatrix.setPrecision(10)
      self.checkarray(sdp.getOutput("cost"),DMatrix(3.20626934048e1),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(3.20626923535e1),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([1.551644595,0.6709672545,0.9814916693,1.406569511,0.9421687787]),digits=5)
      
      self.checkarray(sdp.getOutput("dual"),DMatrix(sp,[2.640261206,0.5605636589,0.5605636589,3.717637107,0.7615505416,-1.513524657,1.139370202,-1.513524657,3.008016978,-2.264413045,1.139370202,-2.264413045,1.704633559,0,0]),digits=5)
      self.checkarray(sdp.getOutput("p"),DMatrix(sp,[0,0,0,0,7.119155551,5.024671489,1.916294752,5.024671489,4.414745792,2.506021978,1.916294752,2.506021978,2.048124139,0.3432465654,4.391169489]),digits=5)

  def test_example2_perm(self):
    self.message("Example2_permuted")
    # Originates from http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
    
    A = DMatrix(0,5)
    c = DMatrix([1.1, -10, 6.6 , 19 , 4.1])

    perm = [5,2,1,0,6,3,4]
    permi = lookupvector(perm,len(perm))
    
    G = -blkdiag([DMatrix([[-1.4,-3.2],[-3.2,-28]]),DMatrix([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0]);

    sp = G.sparsity()
    
    flatdata = [[0.5,5.2,5.2,-5.3,7.8,-2.4,6.0,-2.4,4.2,6.5,6.0,6.5,2.1,-4.5,-3.5],
    [1.7,7.0,7.0,-9.3,-1.9,-0.9,-1.3,-0.9,-0.8,-2.1,-1.3,-2.1,4.0,-0.2,-3.7],
 [6.3,-7.5,-7.5,-3.3,0.2,8.8,5.4,8.8,3.4,-0.4,5.4,-0.4,7.5,-3.3,-4.0],
  [-2.4,-2.5,-2.5,-2.9,3.4,-3.2,-4.5,-3.2,3.0,-4.8,-4.5,-4.8,3.6,4.8,9.7],
  [-6.5,-5.4,-5.4,-6.6,6.7,-7.2,-3.6,-7.2,7.3,-3.0,-3.6,-3.0,-1.4,6.1,-1.5]]

    F = -vertcat([DMatrix(sp,data)[perm,perm] for data in flatdata])
    makeSparse(F)
    
    G = G[perm,perm]
    for sdpsolver, sdp_options in sdpsolvers:
      sdp = sdpsolver(sdpStruct(a=A.sparsity(),g=G.sparsity(),f=F.sparsity()))
      sdp.setOption(sdp_options)
      sdp.init()

      sdp.setInput(G,"g")
      sdp.setInput(c,"c")
      sdp.setInput(F,"f")

      sdp.evaluate()
      DMatrix.setPrecision(10)
      self.checkarray(sdp.getOutput("cost"),DMatrix(3.20626934048e1),digits=5)
      self.checkarray(sdp.getOutput("dual_cost"),DMatrix(3.20626923535e1),digits=5)
      self.checkarray(sdp.getOutput("x"),DMatrix([1.551644595,0.6709672545,0.9814916693,1.406569511,0.9421687787]),digits=5)
      
      self.checkarray(sdp.getOutput("dual")[permi,permi],DMatrix(sp,[2.640261206,0.5605636589,0.5605636589,3.717637107,0.7615505416,-1.513524657,1.139370202,-1.513524657,3.008016978,-2.264413045,1.139370202,-2.264413045,1.704633559,0,0]),digits=5)
      self.checkarray(sdp.getOutput("p")[permi,permi],DMatrix(sp,[0,0,0,0,7.119155551,5.024671489,1.916294752,5.024671489,4.414745792,2.506021978,1.916294752,2.506021978,2.048124139,0.3432465654,4.391169489]),digits=5)

    
    
if __name__ == '__main__':
    unittest.main()

