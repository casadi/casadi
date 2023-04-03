#
#     MIT No Attribution
#
#     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
