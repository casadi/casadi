#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
import casadi as ca
import casadi as c
import numpy
from numpy import random, array, linalg, matrix, zeros, ones, ndarray, eye
import unittest
from types import *
from helpers import *
from copy import deepcopy

import sys
import warnings

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

if ca.has_dple("slicot"):
  dplesolvers.append(("slicot",{"linear_solver": "csparse"}))

def randstable(n,margin=0.8,minimal=0):
  r = margin
  A_ = ca.tril(ca.DM(numpy.random.random((n,n))))
  for i in range(n):
    eig = 0
    while abs(eig)<=minimal:
      eig = numpy.random.random()*r*2-r
    A_[i,i] = eig

  Q = scipy.linalg.orth(numpy.random.random((n,n)))
  return ca.mtimes([Q,A_,Q.T])

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
          
          V_ = [v @ v.T for v in [ca.DM(numpy.random.random((n,n))) for i in range(K)]]
          V2_ = [v @ v.T for v in [ca.DM(numpy.random.random((n,n))) for i in range(K)]]
          S = ca.kron(ca.Sparsity.diag(K),ca.Sparsity.dense(n,n))
          solver = ca.dplesol("solver", Solver,{'a':S,'v':ca.repmat(S,1,2)}, options)
          
          inputs = {"a":ca.dcat(A_), "v": ca.horzcat(ca.dcat(V_),ca.dcat(V2_))}
          
          As = ca.MX.sym("A",S)
          Vs = ca.MX.sym("V",S)
          
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = ca.hcat([(i+i.T)/2 for i in isigma(list(ca.diagsplit(Vs,n))) ])
          
          
          AA = ca.dcat([c.kron(i,i) for i in ca.diagsplit(As,n)])

          A_total = ca.DM.eye(n*n*K) - ca.vertcat(*[AA[-n*n:,:],AA[:-n*n,:]])
          
          
          Pf = ca.solve(A_total,ca.vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))
          P = ca.dcat(ca.horzsplit(P,n))
          
          refsol = ca.Function("refsol", {"a": As,"v":Vs,"p":P},ca.dple_in(),ca.dple_out()).map("map","serial",2,["a"],[])

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
          V_ = [v @ v.T for v in [ca.DM(numpy.random.random((n,n))) for i in range(K)]]
        
          inputs = {"a":ca.hcat(A_), "v": ca.hcat(V_)}
          
          As = ca.MX.sym("A",n,n*K)
          Vs = ca.MX.sym("V",n,n*K)
                    
          def sigma(a):
            return a[1:] + [a[0]]
            
          def isigma(a):
            return [a[-1]] + a[:-1]
          
          Vss = ca.hcat([(i+i.T)/2 for i in isigma(list(ca.horzsplit(Vs,n))) ])
          
          
          AA = ca.dcat([c.kron(i,i) for i in ca.horzsplit(As,n)])

          A_total = ca.DM.eye(n*n*K) - ca.vcat([AA[-n*n:,:],AA[:-n*n,:]])
          
          Pf = ca.solve(A_total,ca.vec(Vss),"csparse")
          P = Pf.reshape((n,K*n))

          solver = ca.Function("solver", {"a": As,"v":Vs,"p":ca.hcat(ca.dplesol(ca.horzsplit(As,n),ca.horzsplit(Vs,n),Solver,options))},ca.dple_in(),ca.dple_out())          
          refsol = ca.Function("refsol", {"a": As,"v":Vs,"p":P},ca.dple_in(),ca.dple_out())

          self.checkfunction(solver,refsol,inputs=inputs,failmessage=str(Solver))
    
if __name__ == '__main__':
    unittest.main()
