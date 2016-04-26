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
import numpy

# Let's construct a block diagonal structure
A = diagcat(1,DM([[2,3],[3,4]]),DM([[5,6,7],[6,8,9],[7,9,10]]),11)
print(A)
A.sparsity().spy()

numpy.random.seed(2)

# We randomly permute this nice structure
perm =  list(numpy.random.permutation(list(range(A.size1()))))
AP = A[perm,perm]

print(AP)
AP.sparsity().spy()

# And use scc to recover the blocks
n,p,r = AP.sparsity().scc()

APrestored = AP[p,p]

print(APrestored)
APrestored.sparsity().spy()
print("# blocks: ", n)
print("block boundaries: ", r[:n])
