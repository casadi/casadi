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
#! Linear solvers
#! =================
#!
#! We demonstrate solving a dense system A.x=b by using different linear solvers.
#!
from casadi import *
from numpy import *
import time

n=100
#$ We generate $A \in \mathbf{R}^{n \times n}$, $x \in \mathbf{R}^{n}$ with $n=100$
A=DMatrix([[cos(i*j)-sin(i) for i in range(n)] for j in range(n)])
x=DMatrix([tan(i) for i in range(n)])

#! We generate the b vector:
b=mul(A,x)

#! We demonstrate the LinearSolver API with CSparse:
s = CSparse(A.sparsity())
s.init()

#! Give it the matrix A
s.setInput(A,0)
#! Do the LU factorization
s.prepare()

#! Give it the matrix b
s.setInput(b,1)

#! And we are off to find x...
s.solve()

x_ = s.output()

#! By looking at the residuals between the x we knew in advance and the computed x, we see that the CSparse solver works
print "Sum of residuals = %.2e" % sumAll(fabs(x-x_))

#! Comparison of different linear solvers
#! ======================================
for name, solver in [("LapackLUDense",LapackLUDense),("LapackQRDense",LapackQRDense),("CSparse",CSparse)]:
  s = solver(A.sparsity()) # We create a solver
  s.init()

  s.setInput(A,0) # Give it the matrix A
  
  t0 = time.time()
  for i in range(100):
    s.prepare()        # Do the LU factorization
  pt = (time.time()-t0)/100

  s.setInput(b,1)  # Give it the matrix b

  t0 = time.time()
  for i in range(100):
    s.solve()
  st = (time.time()-t0)/100
  
  x_ = s.output()

  print ""
  print name
  print "=" * 10
  print "Sum of residuals = %.2e" % sumAll(fabs(x-x_))
  print "Preparation time = %0.2f ms" % (pt*1000)
  print "Solve time       = %0.2f ms" % (st*1000)
  assert(sumAll(fabs(x-x_))<1e-9)
  
#! Note that these 
