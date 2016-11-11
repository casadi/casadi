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
import numpy
from numpy import random, array, linalg, matrix, zeros, ones, ndarray, eye
import unittest
from types import *
from helpers import *
from copy import deepcopy

import sys

if sys.version_info >= (3, 0):
  TupleType = tuple


warnings.filterwarnings("ignore",category=DeprecationWarning)

scipy_available = True
try:
	from scipy.sparse import csr_matrix
	import scipy
	import scipy.linalg
except:
	scipy_available = False

dplesolvers = []

if has_dple("slicot"):
  dplesolvers.append(("slicot",{"linear_solver": "csparse"}))

def randstable(n,margin=0.8,minimal=0):
  r = margin
  A_ = tril(DM(numpy.random.random((n,n))))
  for i in range(n):
    eig = 0
    while abs(eig)<=minimal:
      eig = numpy.random.random()*r*2-r
    A_[i,i] = eig

  Q = scipy.linalg.orth(numpy.random.random((n,n)))
  return mtimes([Q,A_,Q.T])

class Controltests(casadiTestCase):

  @skip(not scipy_available)
  @memory_heavy()
  def test_dple_small(self):
    
    for Solver, options in dplesolvers:
      for K in ([3,4] if args.run_slow else [2,3]):
        for n in [2,3,4]:
          print (Solver, options)
          numpy.random.seed(1)
          print (n,K)
          A_ = [randstable(n) for i in range(K)]
          
          V_ = [mtimes(v,v.T) for v in [DM(numpy.random.random((n,n))) for i in range(K)]]
          V2_ = [mtimes(v,v.T) for v in [DM(numpy.random.random((n,n))) for i in range(K)]]
          S = kron(Sparsity.diag(K),Sparsity.dense(n,n))
          solver = dplesol("solver", Solver,{'a':S,'v':repmat(S,1,2)}, options)
          
          inputs = {"a":dcat(A_), "v": horzcat(dcat(V_),dcat(V2_))}
          
          As = MX.sym("A",S)
          Vs = MX.sym("V",S)
          
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = hcat([(i+i.T)/2 for i in isigma(list(diagsplit(Vs,n))) ])
          
          
          AA = dcat([c.kron(i,i) for i in diagsplit(As,n)])

          A_total = DM.eye(n*n*K) - vertcat(*[AA[-n*n:,:],AA[:-n*n,:]])
          
          
          Pf = solve(A_total,vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))
          P = dcat(horzsplit(P,n))
          
          refsol = Function("refsol", {"a": As,"v":Vs,"p":P},dple_in(),dple_out()).map("map","serial",2,["a"],[])

          self.checkfunction(solver,refsol,inputs=inputs,failmessage=str(Solver))
    
  @skip(not scipy_available)
  @memory_heavy()
  def test_dple_alt_small(self):
    
    for Solver, options in dplesolvers:
      for K in ([3,4] if args.run_slow else [2,3]):
        for n in [2,3,4]:
          print (Solver, options)
          numpy.random.seed(1)
          print (n,K)
          A_ = [randstable(n) for i in range(K)]          
          V_ = [mtimes(v,v.T) for v in [DM(numpy.random.random((n,n))) for i in range(K)]]
        
          inputs = {"a":hcat(A_), "v": hcat(V_)}
          
          As = MX.sym("A",n,n*K)
          Vs = MX.sym("V",n,n*K)
                    
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = hcat([(i+i.T)/2 for i in isigma(list(horzsplit(Vs,n))) ])
          
          
          AA = dcat([c.kron(i,i) for i in horzsplit(As,n)])

          A_total = DM.eye(n*n*K) - vcat([AA[-n*n:,:],AA[:-n*n,:]])
          
          Pf = solve(A_total,vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))

          solver = Function("solver", {"a": As,"v":Vs,"p":hcat(dplesol(horzsplit(As,n),horzsplit(Vs,n),Solver,options))},dple_in(),dple_out())          
          refsol = Function("refsol", {"a": As,"v":Vs,"p":P},dple_in(),dple_out())

          self.checkfunction(solver,refsol,inputs=inputs,failmessage=str(Solver))
    
if __name__ == '__main__':
    unittest.main()
