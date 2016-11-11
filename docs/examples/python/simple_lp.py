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
# -*- coding: utf-8 -*-
from casadi import *

# minimize    3x + 4y
# subject to  x + 2y <= 14
#            3x -  y >= 0
#             x -  y <= 2


# Sparsity of the LP linear term
A = Sparsity.dense(3, 2)

# Create solver
solver = conic('solver', 'qpoases', {'a':A})
#solver = conic('solver', 'clp', {'a':A}) # Use clp

g = DM([3,4])
a = DM([[1, 2],[3, -1], [1, -1]])
lba = DM([-inf, 0, -inf])
uba = DM([14, inf, 2])

sol = solver(g=g, a=a, lba=lba, uba=uba)
print(sol)
