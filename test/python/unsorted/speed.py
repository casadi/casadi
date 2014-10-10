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
from scipy.sparse import *
from casadi import *
from random import *
import time
import numpy as N

"""
Test generating some matrices by random accessing, convert between scipy and casadi types.

Conclusions:
CasADi is significantly quicker than scipy for random accessing and inserting new elements
(about 10 times) when compared to the csr_matrix matrix which uses the same storage format
as CasADi. However, csr_matrix is not intended to be used for random access, a more suitable
scipy class is lil_matrix and then when the assembly is finished, convert back to csr_matrix.
This is somewhat quicker than casadi. Since Matrix<> was designed for flexibility rather
than numerical efficiency, the result is satisfactory.

Converting between scipy csr_matrix and casadi matrix works fine.

Summary:
Use scipy for numerical operations on floating points and casadi when working with symbolic
types.
"""

n = 2000
m = 2000

# ~10 % non-zeros
nel = n*m/100

# tests
solvers = ["scipy","casadi"]

# build with linked lists
use_lil = True

# number of repeats
rep = 4

# Generate test matrices
CM = []
SM = []

for sol in solvers:
  print "testing ", sol
  dur_build = []
  dur_convert = []
  dur_operation = []

  for i in range(rep):
    # create empty matrix
    if sol=="casadi":
      M = DMatrix(n,m)
    else:
      if use_lil:
        M = lil_matrix((n,m))
      else:
        M = csr_matrix((n,m))

    # fill with random elements
    t_before = time.time()
    for i in range(nel):
      row = randrange(n)
      col = randrange(m)
      res = random()
      M[row,col] = res

    # need to convert back to csr
    if sol!="casadi" and use_lil:
      M = csr_matrix(M)
      
    t_build = time.time()
    dur_build.append(t_build-t_before)

    # convert casadi matrix to csr
    if sol=="casadi":
      M = csr_matrix((M.data(),M.col(),M.rowind()),(M.size1(),M.size2()),dtype=float)
    
    # Save the matrix
    SM.append(M)

    # convert scipy to casadi
    dim = M.shape
    col = list(int(i) for i in M.indices)
    rowind = list(int(i) for i in M.indptr)
    M = DMatrix(dim[0],dim[1],col,rowind,M.data)
    
    # Save the matrix
    CM.append(M)

    # convert back again
    if sol!="casadi":
      M = csr_matrix((M.data(),M.col(),M.rowind()),(M.size1(),M.size2()),dtype=float)

    t_convert = time.time()
    dur_convert.append(t_convert-t_build)

    # Make a couple of operations
    if len(SM)>1:
      if sol=="casadi":
        M2_C = dot(CM[0],CM[1])
        #for i in range(100):
          #M = trans(M)
        
        #M2 = trans(M)
        #M3 = trans(M2)
        #M4 = M2 + M3
        #M5 = prod(M4,M)
        #M5 = prod(M4,M)
      else:
        M2_S = N.dot(SM[0],SM[1])
        #for i in range(1000):
          #M = M.transpose()
        #M2 = M.transpose()
        #M3 = M2.transpose()
        #M4 = M2 + M3
        #M5 = N.dot(M4,M)
      
    t_operation = time.time()
    dur_operation.append(t_operation-t_convert)
    

  print "durations to build:                      ", dur_build
  print "durations for converting back and forth: ", dur_convert
  print "durations for operations:                ", dur_operation
  print "number of non-zeros for the last matrix: ",
  if sol=="casadi":
    print M.size()
  else:
    print M.getnnz()


M2_C_S = csr_matrix((M2_C.data(),M2_C.col(),M2_C.rowind()),(M2_C.size1(),M2_C.size2()),dtype=float)
M2_DIFF = M2_C_S-M2_S
print "difference is ", repr(M2_DIFF)