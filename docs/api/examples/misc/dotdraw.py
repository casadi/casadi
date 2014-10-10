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
from casadi.tools import *

#! An SX graph
a = SXElement.sym("a")
b = SXElement.sym("b")

c = sin(a**5 + b)

c = c - b/ sqrt(fabs(c))
print c

dotdraw(c)

#! An SX

dotdraw(SX.sym("x",Sparsity.tril(3)))

dotdraw(SX.sym("x",Sparsity.tril(3))**2)

#! An MX graph
x = MX.sym("x",Sparsity.tril(2))
y = MX.sym("y",Sparsity.tril(2))

z = MX.sym("z",4,2)

zz = x+y

dotdraw(zz)

f = MXFunction([z,y],[z+x[0,0],x-y])
f.setOption("name","magic")
f.init()

[z,z2] = f.call([vertcat([x,y]),zz.T])

z = z[:2,:] +x + cos(x) - sin(x) / tan(z2)

dotdraw(z)
