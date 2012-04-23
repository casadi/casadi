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
#row = [0,0,1,1,2]
#col = [0,0,1,1,2]
#row = [2,1,1,0,0]
#col = [2,0,0,0,0]
row = [0,1,2]
col = [0,0,0]


mapping = IVector()
print "row = ", row
print "col = ", col

sp = sp_triplet(60,1,row,col,mapping,False,True)
print "mapping = ", mapping
print "sp = ", sp
print "sp.col() = ", sp.col()
print "sp.getRow() = ", sp.getRow()

a = DMatrix(10,1)
a[1,0] = 1
nz = IVector([0,1,2])
print "nz = ", nz
a.sparsity().getNZInplace(nz)
print "nz = ", nz
