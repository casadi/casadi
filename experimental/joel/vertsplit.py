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
n = 12
m = 3

x = msym("x",n,m)
row_offset = range(0,n,n/3)
print "row_offset = ", row_offset

r1,r2,r3 = vertsplit(x,row_offset)

f = MXFunction([x],[r3,r2])
f.init()
f.setInput(range(n*m))

d = (n*m-1)*[0] + [1]
f.setFwdSeed(d)

f.evaluate(1,0)
for i in range(2): print f.output(i)
for i in range(2): print f.fwdSens(i)

print f

gf = f.derivative(0,1)
print gf



