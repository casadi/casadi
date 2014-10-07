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
#A = DMatrix.ones(3,3) + DMatrix.ones(sp_diag(3)) # dense
#A = DMatrix.ones(sp_diag(3)) # diagonal
A = 2*DMatrix.ones(sp_diag(3)); A[1,0] = 1; A[2,0] = 1; A[2,1] = 1 # lower triangular

# Permute
A = vertcat((A[1,:],A[2,:],A[0,:]))

print "A = "
A.printDense()

print "SX"
AA = ssym("A",A.sparsity())

r = ssym("r",3)
x = solve(AA,r)
f = SXFunction([r,AA],[x])
f.init()
DMatrix.ones(f.jacSparsity(0,0)).printDense()
DMatrix.ones(f.jacSparsity(1,0)).printDense()

print "SX, implicit function"
x = ssym("x",3)
res = SXFunction([x,r,AA],[mul(AA,x)-r])
f = NewtonImplicitSolver(res)
f.setOption("linear_solver",CSparse)
f.setOption("ad_mode","reverse")
f.init()
DMatrix.ones(f.jacSparsity(0,0)).printDense()
DMatrix.ones(f.jacSparsity(1,0)).printDense()

print "MX"
AA = msym("A",A.sparsity())

r = msym("r",3)
sol = CSparse(AA.sparsity())
sol.init()
x = sol.solve(AA,r.T,True).T
f = MXFunction([r,AA],[x])
#f.setOption("ad_mode","reverse")
f.init()
DMatrix.ones(f.jacSparsity(0,0)).printDense()
DMatrix.ones(f.jacSparsity(1,0)).printDense()

print "MX, implicit function"
x = msym("x",3)
res = MXFunction([x,r,AA],[mul(AA,x)-r])
f = NewtonImplicitSolver(res)
f.setOption("linear_solver",CSparse)
f.setOption("ad_mode","reverse")
f.init()
DMatrix.ones(f.jacSparsity(0,0)).printDense()
DMatrix.ones(f.jacSparsity(1,0)).printDense()


