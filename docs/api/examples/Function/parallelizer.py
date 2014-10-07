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
#! Demonstration of Parallelizer
#! =============================
from casadi import *

n = 12

#! Construct a functon that is expensive to evaluate
x = MX.sym("x",100,100)
y = MX.sym("y",100,1)
z = x

for i in range(20):
  z = mul(x,z)


f = MXFunction([x,y],[z])
f.init()

#! Evaluate this function ten times in parallel
p = Parallelizer([f]*n)
p.setOption("gather_stats",True)
p.setOption("parallelization","openmp")
p.init()

#! Note that the parallelizer MX input/output interface is a repitition of our function's I/O interface
assert(p.getNumInputs() == n*f.getNumInputs())

p.evaluate()

print p.getStats()

