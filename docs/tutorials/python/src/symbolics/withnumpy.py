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
#! CasADi tutorial
#! ==================
#! This tutorial file explains the use of CasADi's SX in a python context.
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *

X=SX.sym('X',2,2)
Y=SX.sym('Y',2,2)
f=SXFunction ([X], [X])
f.init()
print f([Y])

a=SXElement.sym("a")
b=SXElement.sym("b")
c=SXElement.sym("c")
d=SXElement.sym("d")

A=array([[a,b],[c,d]])
B= SX(A)


