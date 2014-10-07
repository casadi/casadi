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
#! CasADi
#! ==================
#! More about SX
from casadi import *
from numpy import *
#! The identity of an SXElement node is very persistant.
#! We demonstrate this with the help of symbolic substitution.
x=SX.sym("x")
y=x**2
f = SXFunction([x],[y])
f.init()
print f([SX.sym("w")])
#! We expect w^2.
l = x
f = SXFunction([l],[y])
f.init()
print f([SX.sym("w")])
#! We expect w^2.
k=SX(x)
l=k[0]
f = SXFunction([l],[y])
f.init()
print f([SX.sym("w")])
#! We expect w^2.
k=SX.sym("d",2,2)
k.nz[1] = x
l=k.nz[1]
f = SXFunction([l],[y])
f.init()
print f([SX.sym("w")])
#! We expect w^2.
#! Identity is not associated with name:
l=SX.sym("x")
f = SXFunction([l],[y])
f.init()
print f([SX.sym("w")])
