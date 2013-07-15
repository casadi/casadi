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

lsolvers = []
try:
  lsolvers.append((CSparse,{}))
except:
  pass
  
try:
  lsolvers.append((LapackLUDense,{}))
except:
  pass
  
try:
  lsolvers.append((LapackQRDense,{}))
except:
  pass
  
#try:
#  lsolvers.append((SymbolicQR,{}))
#except:
#  pass
  
print lsolvers

class LinearSolverTests(casadiTestCase):
  
  def test_simple_trans(self):
    A = DMatrix([[3,7],[1,2]])
    for Solver, options in lsolvers:
      solver = Solver(A.sparsity())
      solver.setOption(options)
      solver.init()
      solver.setInput(A,"A")

      solver.setInput(True,"T")
      solver.prepare()
      
      b = DMatrix([1,0.5])
      solver.setInput(b.T,"B")
      
      solver.solve()
      
      res = DMatrix([1.5,-0.5])
      self.checkarray(solver.output("X"),res.T)
      #   result' = A\b'               Ax = b

  def test_simple(self):
    A = DMatrix([[3,7],[1,2]])
    for Solver, options in lsolvers:
      print Solver
      solver = Solver(A.sparsity())
      solver.setOption(options)
      solver.init()
      solver.setInput(A,"A")

      solver.setInput(False,"T")
      solver.prepare()
      
      b = DMatrix([1,0.5])
      solver.setInput(b.T,"B")
      
      solver.solve()
      
      res = DMatrix([-1.5,5.5])
      self.checkarray(solver.output("X"),res.T)
      #   result' = A'\b'             Ax = b
      
      
  @requires("CSparseCholesky")
  def test_cholesky(self):
    random.seed(1)
    n = 10
    L = self.randDMatrix(n,n,sparsity=0.2) +  c.diag(range(n))
    M = mul(L,L.T)

    S = CSparseCholesky(M.sparsity())

    S.init()
    S.getFactorizationSparsity().spy()

    S.setInput(M,0)

  @requires("CSparseCholesky")
  def test_cholesky2(self):
    random.seed(0)
    n = 10
    L = c.diag(range(n))
    M = mul(L,L.T)

    print L
    S = CSparseCholesky(M.sparsity())
    

    S.init()
    S.getFactorizationSparsity().spy()

if __name__ == '__main__':
    unittest.main()
