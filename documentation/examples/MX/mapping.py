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
#! createParent
#!======================
from casadi import *
from numpy import *
import casadi as c

##! We create a bunch of matrices that depend on V
V = MX("V",13)
A = c.reshape(V[:2],sp_dense(2,1))  # 2x1 matrix
B = c.reshape(V[2:3],sp_dense(1,1)) # 1x1 matrix
C = c.reshape(V[3:],sp_tril(4))     # Triangular matrix

print "A = ", A
print "A.mapping() = ", A.mapping()
print "B = ", B
print "B.mapping() = ", B.mapping()
print "C = ", C
print "C.mapping() = ", C.mapping()

##! We create a bunch of matrices that depend on V
V_= DMatrix(13,1,True)   # create a dense DMatrix

##! We can use mapping() to index into this DMatrix:
V_[A.mapping()] = 1
V_[B.mapping()] = 10
V_[C.mapping()] = 100
  
print V_
