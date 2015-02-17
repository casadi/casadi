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
import random
import time

import __builtin__
import scipy.linalg

clesolvers = []
if LinearSolver.hasPlugin("csparse") and CleSolver.hasPlugin("simple"):
  clesolvers.append(("simple",{"linear_solver": "csparse"}))

dlesolvers = []

if LinearSolver.hasPlugin("csparse") and DleSolver.hasPlugin("simple"):
  dlesolvers.append(("simple",{"linear_solver": "csparse"}))
  
  
if LinearSolver.hasPlugin("csparse") and DleSolver.hasPlugin("dple.slicot"):
  dlesolvers.append(("dple.slicot",{"dple_solver_options": {"linear_solver": "csparse"}}))

"""
if DleSolver.hasPlugin("lrdle.smith"):
  dlesolvers.append(("lrdle.smith",{"lrdle_solver_options": {"max_iter":100,"tol": 1e-13}}))
"""

if DleSolver.hasPlugin("lrdle.fixed_smith"):
  dlesolvers.append(("lrdle.fixed_smith",{"lrdle_solver_options": {"iter":100}}))

  
if DleSolver.hasPlugin("fixed_smith"):
  dlesolvers.append(("fixed_smith",{"iter":100, "freq_doubling": False}))
  dlesolvers.append(("fixed_smith",{"iter":100, "freq_doubling": True}))

lrdlesolvers = []

"""
if LrDleSolver.hasPlugin("smith"):
  lrdlesolvers.append(("smith",{"max_iter":100,"tol": 1e-13}))
"""

if LrDleSolver.hasPlugin("fixed_smith"):
  lrdlesolvers.append(("fixed_smith",{"iter":100}))


if LrDleSolver.hasPlugin("dle.simple") and LinearSolver.hasPlugin("csparse"):
  lrdlesolvers.append(("dle.simple",{"dle_solver_options": {"linear_solver": "csparse"}}))

"""
if LrDleSolver.hasPlugin("dle.lrdle.smith"):
  lrdlesolvers.append(("dle.lrdle.smith",{"dle_solver_options.lrdle_solver_options": {"max_iter":100,"tol": 1e-13}}))
"""

lrdplesolvers = []

""" lrsmith has a problem #1294
if LrDpleSolver.hasPlugin("lifting.smith"):
  lrdplesolvers.append(("lifting.smith",{"lrdle_solver_options": {"max_iter":100,"tol": 1e-13}}))
  lrdplesolvers.append(("lifting.smith",{"form": "B", "lrdle_solver_options": {"max_iter":100,"tol": 1e-13}}))
"""

if LrDpleSolver.hasPlugin("dple.slicot") and LinearSolver.hasPlugin("csparse"):
  lrdplesolvers.append(("dple.slicot",{"dple_solver_options": {"linear_solver": "csparse"}}))

""" lrsmith has a problem #1294
if LrDpleSolver.hasPlugin("dple.lrdple.lifting.smith"):
  lrdplesolvers.append(("dple.lrdple.lifting.smith",{ "dple_solver_options.lrdple_solver_options.lrdle_solver_options": {"max_iter":100,"tol": 1e-13}} ))
"""

dplesolvers = []

if LinearSolver.hasPlugin("csparse") and DpleSolver.hasPlugin("slicot"):
  dplesolvers.append(("slicot",{"linear_solver": "csparse"}))
  
if LinearSolver.hasPlugin("csparse") and DpleSolver.hasPlugin("simple"):
  dplesolvers.append(("simple",{"linear_solver": "csparse"}))

if LinearSolver.hasPlugin("csparse") and DpleSolver.hasPlugin("condensing.simple"):
  dplesolvers.append(("condensing.simple",{"dle_solver_options": {"linear_solver": "csparse"}}))

if DpleSolver.hasPlugin("condensing.dple.slicot") and LinearSolver.hasPlugin("csparse"):
  dplesolvers.append(("condensing.dple.slicot",{"dle_solver_options.dple_solver_options": {"linear_solver": "csparse"}}))

if DpleSolver.hasPlugin("lrdple.lifting.dle.dple.slicot") and LinearSolver.hasPlugin("csparse"):
  dplesolvers.append(("lrdple.lifting.dle.dple.slicot",{"lrdple_solver_options.lrdle_solver_options.dle_solver_options.dple_solver_options": {"linear_solver": "csparse"}}))

if DpleSolver.hasPlugin("lifting.dple.simple") and LinearSolver.hasPlugin("csparse"):
  dplesolvers.append(("lifting.dple.simple",{"dle_solver_options.dple_solver_options": {"linear_solver": "csparse"}}))


if DpleSolver.hasPlugin("lifting.dple.slicot") and LinearSolver.hasPlugin("csparse"):
  dplesolvers.append(("lifting.dple.slicot",{"dle_solver_options.dple_solver_options": {"linear_solver": "csparse"}}))
  dplesolvers.append(("lifting.dple.slicot",{"form": "B", "dle_solver_options.dple_solver_options": {"linear_solver": "csparse"}}))

print "DpleSolvers", len(dplesolvers), dplesolvers
print "DleSolvers", len(dlesolvers), dlesolvers
print "CleSolvers", len(clesolvers), clesolvers


print "LrDpleSolvers", len(lrdplesolvers),  lrdplesolvers
print "LrDleSolvers", len(lrdlesolvers), lrdlesolvers

if not args.run_slow:
  dplesolvers = dplesolvers[:1]
  dlesolvers = dlesolvers[:1]
  clesolvers = clesolvers[:1]
  lrdplesolvers = lrdplesolvers[:1]
  lrdlesolvers = lrdlesolvers[:1]
  
def randstable(n,margin=0.8):
  r = margin
  A_ = tril(DMatrix(numpy.random.random((n,n))))
  for i in range(n): A_[i,i] = numpy.random.random()*r*2-r

  Q = scipy.linalg.orth(numpy.random.random((n,n)))
  return mul([Q,A_,Q.T])

class ControlTests(casadiTestCase):
  
  @slow()
  @memory_heavy()
  def test_dple_CH(self):

    for with_C in [True]:
      for with_H in [True]:
        print "with_C", with_C
        print "with_H", with_H
        
        numpy.random.seed(1)
        
        n = 5
        
        m = 3 if with_C else n
        h = [3,1]
        
        Ls = [0]+list(numpy.cumsum(h))
         
        K = 3

        H_ = [DMatrix(numpy.random.random((n,__builtin__.sum(h))))  if with_H else DMatrix() for i in range(K)]

        A_ = [ randstable(n,margin=0.2) for i in range(K) ]
  
        C_ = [ DMatrix(numpy.random.random((n,m))) if with_C else DMatrix() for i in range(K) ]
        V_ = [ DMatrix(numpy.random.random((m,m))) for i in range(K) ]
        V_ = [ (i+i.T)/2 for i in V_ ] 
        
        H = MX.sym("H",horzcat(H_).sparsity())
        Hs_ = horzsplit(H,__builtin__.sum(h))
        Hss_ = [ horzsplit(i,Ls) for i in Hs_]
        
        A = MX.sym("A",horzcat(A_).sparsity())
        As_ = horzsplit(A,n)
        
        C = MX.sym("C",horzcat(C_).sparsity())
        Cs_ = horzsplit(C,m)
        
        Vs = V = MX.sym("V",horzcat(V_).sparsity())
        Vs_ = [(i+i.T)/2 for i in horzsplit(V,m)]
        
        print "Vs_", [i.dimString() for i in Vs_]
        #V = horzcat([ (i+i.T)/2 for i in Vs_])

        N = 100
        
        V = diagcat([Vs_[-1],diagcat(Vs_[:-1])])
        
        Hn = diagcat(Hs_)
        An = blockcat([[MX(n,(K-1)*n),As_[-1]],[diagcat(As_[:-1]),MX((K-1)*n,n)]])
        Cn = diagcat([Cs_[-1],diagcat(Cs_[:-1])]) if with_C else None
        Vn = diagcat([Vs_[-1],diagcat(Vs_[:-1])])

        D = [Cn if with_C else DMatrix.eye(n*K)]
        for i in range(N):
          D.append(mul(An,D[-1]))

        DD = horzcat(D)
        
        Ls2 = [i*K for i in Ls]
        
        Ls2 = [0]+list(numpy.cumsum(h*K))

      
        if with_H:
          Y = [ diagcat([ mul([mul(Li.T,DD),diagcat([Vn]*(N+1)),mul(DD.T,Li)]) for Li in horzsplit(Hnn,Ls)]) for Hnn in horzsplit(Hn,__builtin__.sum(h))]
        else: 
          Y = diagsplit(mul([DD,diagcat([Vn]*(N+1)),DD.T]),n)
        
        
        if with_C:
          if with_H:
            for i in Y:
              i.sparsity().spy()
            f = MXFunction(lrdpleIn(a=A,c=C,v=Vs,h=H),[horzcat(Y)])
          else:
            f = MXFunction(lrdpleIn(a=A,c=C,v=Vs),[horzcat(Y)])
        else:
          if with_H:
            f = MXFunction(lrdpleIn(a=A,v=Vs,h=H),[horzcat(Y)])
          else:
            f = MXFunction(lrdpleIn(a=A,v=Vs),[horzcat(Y)])
        f.setOption("name","reference")
        f.init()
        
        print "f",f
        
        
        temp = MXFunction([A],[An])
        temp.init()
        Af_ = temp([horzcat(A_)])[0]
        print "Af"
        Af_.printDense()
        E = numpy.linalg.eig(Af_)[0]
        
        assert max(abs(E))<=0.95, str(max(abs(E)))
        
        
        
        print [ numpy.linalg.eig(i)[0] for i in A_ ] 
        
        print numpy.linalg.eig(mul([i.T for i in reversed(A_)]))[0]
        print numpy.linalg.eig(mul([i for i in reversed(A_)]))[0]
        print numpy.linalg.eig(mul([i.T for i in A_]))[0]
        print numpy.linalg.eig(mul([i for i in A_]))[0]
        
        
        for Solver, options in lrdplesolvers:
          print Solver, options
          print "c", [i.sparsity() for i in Cs_]
          g = LrDpleSolver(Solver,lrdpleStruct(a=[i.sparsity() for i in As_],c=[i.sparsity() for i in Cs_],v=[ i.sparsity() for i in Vs_],h=[i.sparsity() for i in Hs_]if with_H else [])  )
          
          print g.dictionary()
          g.setOption("Hs",[h]*K)
          print options
          g.setOption(options)
          g.init()
          
          for i in [f,g]:
            i.setInput(horzcat(A_),"a")
            i.setInput(horzcat(V_),"v")
            if with_C:
              i.setInput(horzcat(C_),"c")
            if with_H:
              i.setInput(horzcat(H_),"h")
              
            i.evaluate()
                 
                    
          try:
            self.checkfunction(g,f,sens_der=True,hessian=True,evals=2)
          except Exception as e:
            if "second order derivatives are not supported" in str(e):
              self.checkfunction(g,f,evals=1,hessian=False,sens_der=False)
            else:
              raise e

  @memory_heavy()
  def test_custom2(self):
    numpy.random.seed(1)
    
    for with_C in [True]:
      for with_H in [True]:
        print "with_C", with_C
        print "with_H", with_H
        
        n = 5
        
        m = 3 if with_C else n
        h = [2,3]
        h = [3,1]

        H_ = DMatrix(numpy.random.random((n,__builtin__.sum(h))))  if with_H else DMatrix()
        P_ = numpy.random.random((n,n))
        R_ = numpy.random.random((n,3))

        A_ = randstable(n)

        C_ = DMatrix(numpy.random.random((n,m))) if with_C else DMatrix()
        V_ = DMatrix(numpy.random.random((m,m)))
        V_ = (V_+V_.T)/2
        
        H = MX.sym("H",H_.sparsity())

        A = MX.sym("A",A_.sparsity())
        C = MX.sym("C",C_.sparsity())

        Vs = V = MX.sym("V",V_.sparsity())
        V = (V+V.T)/2

        N = 100

        D = [C if with_C else DMatrix.eye(n)]
        for i in range(N):
          D.append(mul(A,D[-1]))

        DD = horzcat(D)
        Ls = [0]+list(numpy.cumsum(h))
      
        if with_H:
          Y = [ mul([mul(Li.T,DD),diagcat([V]*(N+1)),mul(DD.T,Li)]) for Li in horzsplit(H,Ls)]
        else: 
          Y = [ mul([DD,diagcat([V]*(N+1)),DD.T]) ]
        
        if with_C:
          if with_H:
            f = MXFunction(lrdleIn(a=A,c=C,v=Vs,h=H),[diagcat(Y)])
          else:
            f = MXFunction(lrdleIn(a=A,c=C,v=Vs),[diagcat(Y)])
        else:
          if with_H:
            f = MXFunction(lrdleIn(a=A,v=Vs,h=H),[diagcat(Y)])
          else:
            f = MXFunction(lrdleIn(a=A,v=Vs),[diagcat(Y)])
        f.init()
        
        
        for Solver, options in lrdlesolvers:
          print Solver
          g = LrDleSolver(Solver,lrdleStruct(a=A.sparsity(),c=C.sparsity(),v=Vs.sparsity(),h=H.sparsity()))
          g.setOption("Hs",h  if with_H else [])
          g.setOption(options)
          g.init()
          
          for i in [f,g]:
            i.setInput(A_,"a")
            i.setInput(V_,"v")
            if with_C:
              i.setInput(C_,"c")
            if with_H:
              i.setInput(H_,"h")
          
          try:
            self.checkfunction(g,f,sens_der=True,hessian=True,evals=2,digits=7)
          except Exception as e:
            if "second order derivatives are not supported" in str(e):
              self.checkfunction(g,f,evals=1,hessian=False,sens_der=False,digits=7)
            else:
              raise e
  @memory_heavy()
  def test_dle_small(self):
    
    for Solver, options in dlesolvers:
        for n in [2,3,4]:
          numpy.random.seed(1)
          print (n)
          
          r = 0.8
          A_ = tril(DMatrix(numpy.random.random((n,n))))
          for i in range(n): A_[i,i] = numpy.random.random()*r*2-r
          
          Q = scipy.linalg.orth(numpy.random.random((n,n)))
          A_ = mul([Q,A_,Q.T])
          
          v = DMatrix(numpy.random.random((n,n)))
          V_ = mul(v,v.T)
          
          
          solver = DleSolver(Solver,dleStruct(a=Sparsity.dense(n,n),v=Sparsity.dense(n,n)))
          solver.setOption(options)
          solver.init()
          solver.setInput(A_,"a")
          solver.setInput(V_,"v")
          
          As = MX.sym("A",n,n)
          Vs = MX.sym("V",n,n)
          
          Vss = (Vs+Vs.T)/2
          
          A_total = DMatrix.eye(n*n) - c.kron(As,As)
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          
          refsol = MXFunction(dleIn(a=As,v=Vs),dleOut(p=Pf.reshape((n,n))))
          refsol.init()
          
          refsol.setInput(A_,"a")
          refsol.setInput(V_,"v")
          
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
          
          try:
            self.checkfunction(solver,refsol,sens_der=True,hessian=True,evals=2,failmessage=str(Solver))
          except Exception as e:
            if "second order derivatives are not supported" in str(e):
              self.checkfunction(solver,refsol,evals=1,hessian=False,sens_der=False,failmessage=str(Solver))
            else:
              raise e
       
  @memory_heavy()
  def test_cle_small(self):
    
    for Solver, options in clesolvers:
        for n in [2,3,4]:
          numpy.random.seed(1)
          print (n)
          A_ = DMatrix(numpy.random.random((n,n)))
          
          v = DMatrix(numpy.random.random((n,n)))
          V_ = mul(v,v.T)
          
          
          solver = CleSolver(Solver,cleStruct(a=Sparsity.dense(n,n),v=Sparsity.dense(n,n)))
          solver.setOption(options)
          solver.init()
          solver.setInput(A_,"a")
          solver.setInput(V_,"v")
          
          As = MX.sym("A",n,n)
          Vs = MX.sym("V",n,n)
          
          Vss = (Vs+Vs.T)/2
          
          e = DMatrix.eye(n)
          
          A_total = - c.kron(e,As) - c.kron(As,e)
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          
          refsol = MXFunction(dleIn(a=As,v=Vs),dleOut(p=Pf.reshape((n,n))))
          refsol.init()
          
          refsol.setInput(A_,"a")
          refsol.setInput(V_,"v")
          
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
          print Solver, options
          numpy.random.seed(1)
          print (n,K)
          A_ = [randstable(n) for i in range(K)]
          
          V_ = [mul(v,v.T) for v in [DMatrix(numpy.random.random((n,n))) for i in range(K)]]
          
          
          solver = DpleSolver(Solver,dpleStruct(a=[Sparsity.dense(n,n) for i in range(K)],v=[Sparsity.dense(n,n) for i in range(K)]))
          solver.setOption(options)
          solver.init()
          solver.setInput(horzcat(A_),"a")
          solver.setInput(horzcat(V_),"v")
          
          As = MX.sym("A",n,K*n)
          Vs = MX.sym("V",n,K*n)
          
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = horzcat([(i+i.T)/2 for i in isigma(list(horzsplit(Vs,n))) ])
          
          
          AA = diagcat([c.kron(i,i) for i in horzsplit(As,n)])

          A_total = DMatrix.eye(n*n*K) - vertcat([AA[-n*n:,:],AA[:-n*n,:]])
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))
          
          refsol = MXFunction(dpleIn(a=As,v=Vs),dpleOut(p=P))
          refsol.init()
          
          refsol.setInput(horzcat(A_),"a")
          refsol.setInput(horzcat(V_),"v")
          
          solver.evaluate()
          X = list(horzsplit(solver.getOutput(),n))
          refsol.evaluate()
          Xref = list(horzsplit(refsol.getOutput(),n))
          
          a0 = (mul([diagcat(A_),diagcat(X),diagcat(A_).T])+diagcat(V_))
          a0ref = (mul([diagcat(A_),diagcat(Xref),diagcat(A_).T])+diagcat(V_))
          

            
          a1 = diagcat(sigma(X))
          a1ref = diagcat(sigma(Xref))

          self.checkarray(a0ref,a1ref,failmessage=str(Solver))
          self.checkarray(a0,a1,failmessage=str(Solver))
          
          try:
            self.checkfunction(solver,refsol,sens_der=True,hessian=True,evals=2,failmessage=str(Solver))
          except Exception as e:
            if "second order derivatives are not supported" in str(e):
              self.checkfunction(solver,refsol,evals=1,hessian=False,sens_der=False,failmessage=str(Solver))
            else:
              raise e
  
  @memory_heavy()
  def test_dple_large(self):
    
    for Solver, options in dplesolvers:  
      for K in ([1,2,3,4,5] if args.run_slow else [1,2,3]):
        for n in ([2,3,4,8,16,32] if args.run_slow else [2,3,4]):
          if ("simple" in str(Solver) or "simple" in str(options)) and n*K > 40 : continue
          if "lifting" in Solver and "dple.slicot" in Solver+str(options)  and K>4: continue
          # Somehow, for K>4 we are constructing a system that slicot fails to solve properly
          
          print Solver, options
          numpy.random.seed(1)
          print (n,K)
          A_ = [randstable(n) for i in range(K)]
          
          V_ = [mul(v,v.T) for v in [DMatrix(numpy.random.random((n,n))) for i in range(K)]]
          
          
          solver = DpleSolver(Solver,dpleStruct(a=[Sparsity.dense(n,n) for i in range(K)],v=[Sparsity.dense(n,n) for i in range(K)]))
          solver.setOption(options)
          solver.init()
          solver.setInput(horzcat(A_),"a")
          solver.setInput(horzcat(V_),"v")
          
          t0 = time.time()
          solver.evaluate()
          print "eval [ms]: ", (time.time()-t0)*1000
          X = list(horzsplit(solver.getOutput(),n))

          def sigma(a):
            return a[1:] + [a[0]]
            
          for a,v,x,xp in zip(A_,V_,X,sigma(X)):
            self.checkarray(xp,mul([a,x,a.T])+v,digits=2 if "condensing" in str(Solver) else 7)
          
  @requiresPlugin(DpleSolver,"slicot")
  def test_slicot_periodic_schur(self):
    for K in ([1,2,3,4,5] if args.run_slow and not args.ignore_memory_heavy else [1,2,3]):
      for n in ([2,3,4,8,16,32] if args.run_slow and not args.ignore_memory_heavy else [2,3,4]):
        numpy.random.seed(1)
        A = [DMatrix(numpy.random.random((n,n))) for i in range(K)]
        T,Z,er,ec = DpleSolver.periodic_schur('slicot',A)
        def sigma(a):
          return a[1:] + [a[0]]
              
        for z,zp,a,t in zip(Z,sigma(Z),A,T):
          self.checkarray(mul([z.T,a,zp]),t,digits=7)
          
        hess = Sparsity.band(n,1)+Sparsity.upper(n)
        # T[0]  in hessenberg form
        self.checkarray(T[0][hess.patternInverse()],DMatrix.zeros(n,n),digits=12)
        
        # remainder of T is upper triangular
        for t in T[1:]:
          self.checkarray(t[Sparsity.upper(n).patternInverse()],DMatrix.zeros(n,n),digits=12)
          
        for z in Z:
          self.checkarray(mul(z,z.T),DMatrix.eye(n))
          self.checkarray(mul(z.T,z),DMatrix.eye(n))
          
if __name__ == '__main__':
    print sys.argv
    unittest.main()
