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

# This pattern should work fine for unidirectional coloring
A = DMatrix.eye(5)
A[0,:] = 1

# Unidirectional coloring
coloring = A.sparsity().unidirectionalColoring()
print coloring

# Create a symmetric matrix with a for unidirectional coloring "bad" sparsity pattern
A = DMatrix.eye(5)
#A = DMatrix.zeros(5,5)
A[:,-1] = 1
A[-1,:] = 1

# Unidirectional coloring
coloring = A.sparsity().unidirectionalColoring()
print coloring

# Star coloring
coloring = A.sparsity().starColoring()
print coloring

# Largest first ordering
ordering = A.sparsity().largestFirstOrdering()
print "ordering = ", ordering

# Create a function whose hessian has the corresponding sparsity pattern
x = ssym("x",5)
y = sin(x[4])
f = y*(x[0]+x[1]+x[2]+x[3])
ff = SXFunction([x],[f])
ff.init()

gf = ff.jac().trans()
#gf = gf[0]

gff = SXFunction([x],[gf])
gff.init()

hf = gff.jac(0,0,False,True)
hf.printDense()

hff2 = ff.hess()
hff2.printDense()



# Unidi

