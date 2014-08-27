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
import random
import time

clesolvers = []
try:
  LinearSolver.loadPlugin("csparse")
  CleSolver.loadPlugin("simple")
  clesolvers.append(("simple",{"linear_solver": "csparse"}))
except:
  pass


dlesolvers = []
try:
  LinearSolver.loadPlugin("csparse")
  DleSolver.loadPlugin("simple")
  dlesolvers.append(("simple",{"linear_solver": "csparse"}))
except:
  pass

try:
  LinearSolver.loadPlugin("csparse")
  DleSolver.loadPlugin("dple")
  dlesolvers.append(("dple",{"dple_solver": "simple", "dple_solver_options": {"linear_solver": "csparse"}}))
except:
  pass
  
dplesolvers = []
try:
  LinearSolver.loadPlugin("csparse")
  DpleSolver.loadPlugin("slicot")
  dplesolvers.append(("slicot",{"linear_solver": "csparse"}))
except:
  pass
  
  
try:
  LinearSolver.loadPlugin("csparse")
  DpleSolver.loadPlugin("simple")
  dplesolvers.append(("simple",{"linear_solver": "csparse"}))
except:
  pass
  
try:
  LinearSolver.loadPlugin("csparse")
  DpleSolver.loadPlugin("condensing")
  dplesolvers.append(("condensing",{"dle_solver": "simple","dle_solver_options": {"linear_solver": "csparse"}}))
except:
  pass
  
print "DpleSolvers", dplesolvers
print "DleSolvers", dlesolvers
print "CleSolvers", clesolvers

class ControlTests(casadiTestCase):

  @memory_heavy()
  def test_dle_small(self):
    
    for Solver, options in dlesolvers:
        for n in [2,3,4]:
          numpy.random.seed(1)
          print (n)
          A_ = DMatrix(numpy.random.random((n,n)))
          
          v = DMatrix(numpy.random.random((n,n)))
          V_ = mul(v,v.T)
          
          
          solver = DleSolver(Solver,Sparsity.dense(n,n),Sparsity.dense(n,n))
          solver.setOption(options)
          solver.init()
          solver.setInput(A_,DLE_A)
          solver.setInput(V_,DLE_V)
          
          As = MX.sym("A",n,n)
          Vs = MX.sym("V",n,n)
          
          Vss = (Vs+Vs.T)/2
          
          A_total = DMatrix.eye(n*n) - c.kron(As,As)
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          
          refsol = MXFunction([As,Vs],[Pf.reshape((n,n))])
          refsol.init()
          
          refsol.setInput(A_,DLE_A)
          refsol.setInput(V_,DLE_V)
          
          solver.evaluate()
          X = solver.getOutput()
          refsol.evaluate()
          Xref = refsol.getOutput()
          
          a0 = (mul([A_,X,A_.T])+V_)
          a0ref = (mul([A_,Xref,A_.T])+V_)
          

            
          a1 = X
          a1ref = Xref
          
          self.checkarray(a0ref,a1ref)
          self.checkarray(a0,a1)

          self.checkfunction(solver,refsol,sens_der=True,hessian=True,evals=2)

  @memory_heavy()
  def test_cle_small(self):
    
    for Solver, options in clesolvers:
        for n in [2,3,4]:
          numpy.random.seed(1)
          print (n)
          A_ = DMatrix(numpy.random.random((n,n)))
          
          v = DMatrix(numpy.random.random((n,n)))
          V_ = mul(v,v.T)
          
          
          solver = CleSolver(Solver,Sparsity.dense(n,n),Sparsity.dense(n,n))
          solver.setOption(options)
          solver.init()
          solver.setInput(A_,DLE_A)
          solver.setInput(V_,DLE_V)
          
          As = MX.sym("A",n,n)
          Vs = MX.sym("V",n,n)
          
          Vss = (Vs+Vs.T)/2
          
          e = DMatrix.eye(n)
          
          A_total = - c.kron(e,As) - c.kron(As,e)
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          
          refsol = MXFunction([As,Vs],[Pf.reshape((n,n))])
          refsol.init()
          
          refsol.setInput(A_,DLE_A)
          refsol.setInput(V_,DLE_V)
          
          solver.evaluate()
          X = solver.getOutput()
          refsol.evaluate()
          Xref = refsol.getOutput()
          
          a0 = mul([A_,X]) + mul([X,A_.T])+V_
          a0ref = mul([A_,Xref]) + mul([Xref,A_.T])+V_
          
          self.checkarray(a0,a0ref)
          self.checkarray(a0,DMatrix.zeros(n,n))

          self.checkfunction(solver,refsol,sens_der=True,hessian=True,evals=2)

  
  @memory_heavy()
  def test_dple_small(self):
    
    for Solver, options in dplesolvers:
      for K in ([1,2,3,4] if args.run_slow else [1,2,3]):
        for n in [2,3]:
          numpy.random.seed(1)
          print (n,K)
          A_ = [DMatrix(numpy.random.random((n,n))) for i in range(K)]
          
          V_ = [mul(v,v.T) for v in [DMatrix(numpy.random.random((n,n))) for i in range(K)]]
          
          
          solver = DpleSolver(Solver,[Sparsity.dense(n,n) for i in range(K)],[Sparsity.dense(n,n) for i in range(K)])
          solver.setOption(options)
          solver.init()
          solver.setInput(horzcat(A_),DPLE_A)
          solver.setInput(horzcat(V_),DPLE_V)
          
          As = MX.sym("A",n,K*n)
          Vs = MX.sym("V",n,K*n)
          
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = horzcat([(i+i.T)/2 for i in isigma(list(horzsplit(Vs,n))) ])
          
          
          AA = blkdiag([c.kron(i,i) for i in horzsplit(As,n)])

          A_total = DMatrix.eye(n*n*K) - vertcat([AA[-n*n:,:],AA[:-n*n,:]])
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))
          
          refsol = MXFunction([As,Vs],[P])
          refsol.init()
          
          refsol.setInput(horzcat(A_),DPLE_A)
          refsol.setInput(horzcat(V_),DPLE_V)
          
          solver.evaluate()
          X = list(horzsplit(solver.getOutput(),n))
          refsol.evaluate()
          Xref = list(horzsplit(refsol.getOutput(),n))
          
          a0 = (mul([blkdiag(A_),blkdiag(X),blkdiag(A_).T])+blkdiag(V_))
          a0ref = (mul([blkdiag(A_),blkdiag(Xref),blkdiag(A_).T])+blkdiag(V_))
          

            
          a1 = blkdiag(sigma(X))
          a1ref = blkdiag(sigma(Xref))

          self.checkarray(a0ref,a1ref)
          self.checkarray(a0,a1)

          self.checkfunction(solver,refsol,sens_der=True,hessian=True,evals=2)
  
  @memory_heavy()
  def test_dple_large(self):
    
    for Solver, options in dplesolvers:
      if "simple" in str(Solver) or "simple" in str(options): continue
      for K in ([1,2,3,4,5] if args.run_slow else [1,2,3]):
        for n in ([2,3,4,8,16,32] if args.run_slow else [2,3,4]):
          numpy.random.seed(1)
          print (n,K)
          A_ = [DMatrix(numpy.random.random((n,n))) for i in range(K)]
          
          V_ = [mul(v,v.T) for v in [DMatrix(numpy.random.random((n,n))) for i in range(K)]]
          
          
          solver = DpleSolver(Solver,[Sparsity.dense(n,n) for i in range(K)],[Sparsity.dense(n,n) for i in range(K)])
          solver.setOption(options)
          solver.init()
          solver.setInput(horzcat(A_),DPLE_A)
          solver.setInput(horzcat(V_),DPLE_V)
          
          t0 = time.time()
          solver.evaluate()
          print "eval [ms]: ", (time.time()-t0)*1000
          X = list(horzsplit(solver.getOutput(),n))

          def sigma(a):
            return a[1:] + [a[0]]
            
          for a,v,x,xp in zip(A_,V_,X,sigma(X)):
            self.checkarray(xp,mul([a,x,a.T])+v,digits=7)
          
  @requires("slicot_periodic_schur")
  def test_slicot_periodic_schur(self):
    for K in ([1,2,3,4,5] if args.run_slow and not args.ignore_memory_heavy else [1,2,3]):
      for n in ([2,3,4,8,16,32] if args.run_slow and not args.ignore_memory_heavy else [2,3,4]):
        numpy.random.seed(1)
        A = [DMatrix(numpy.random.random((n,n))) for i in range(K)]
        T,Z,er,ec = slicot_periodic_schur(A)
        def sigma(a):
          return a[1:] + [a[0]]
              
        for z,zp,a,t in zip(Z,sigma(Z),A,T):
          self.checkarray(mul([z.T,a,zp]),t,digits=7)
          
        hess = Sparsity.band(n,1)+Sparsity.triu(n)
        # T[0]  in hessenberg form
        self.checkarray(T[0][hess.patternInverse()],DMatrix.zeros(n,n),digits=12)
        
        # remainder of T is upper triangular
        for t in T[1:]:
          self.checkarray(t[Sparsity.triu(n).patternInverse()],DMatrix.zeros(n,n),digits=12)
          
        for z in Z:
          self.checkarray(mul(z,z.T),DMatrix.eye(n))
          self.checkarray(mul(z.T,z),DMatrix.eye(n))
          
if __name__ == '__main__':
    unittest.main()
