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

# Let's construct a block diagonal structure
A = blkdiag([1,DMatrix([[2,3],[4,5]]),DMatrix([[6,7,8],[9,10,11],[12,13,14]]),15])
A.printMatrix()

rowperm = IVector()
colperm = IVector()
rowblock = IVector()
colblock = IVector()
coarse_rowblock = IVector()
coarse_colblock = IVector()

A.sparsity().dulmageMendelsohn 	( rowperm, colperm, rowblock,	colblock,	coarse_rowblock, coarse_colblock)

print "rowperm: ", rowperm
print "colperm: ", colperm
print "rowblock: ", rowblock
print "colblock: ", colblock
print "coarse_rowblock: ", coarse_rowblock
print "coarse_colblock: ", coarse_colblock

# Ruin the nice structure

A[0,-1] = 666
A[0,-2] = 666
A.printMatrix()
A.sparsity().dulmageMendelsohn 	( rowperm, colperm, rowblock,	colblock,	coarse_rowblock, coarse_colblock)

print "rowperm: ", rowperm
print "colperm: ", colperm
print "permuted:"
A[rowperm,colperm].printMatrix()
print "rowblock: ", rowblock
print "colblock: ", colblock
print "coarse_rowblock: ", coarse_rowblock
print "coarse_colblock: ", coarse_colblock
