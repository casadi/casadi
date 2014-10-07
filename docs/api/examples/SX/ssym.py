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
#! SX.sym
#!======================
from casadi import *

#! Construct using a single name
#! =====================================
#! The names of the entries of the SX will be derived from the name provided as argument to SX.sym.

#! Without shape arguments, a 1-by-1 matrix is constructed:

x = SX.sym("x")
print type(x), x

#! Create a column matrix 
print SX.sym("x",2,1)

#! Create a row matrix 
print SX.sym("x",1,2)

#! Create a matrix 
print SX.sym("x",2,3)

#! Construct using multiple names
#! =====================================

#! casADi defaults to a column matrix when no shape arguments are provided
print SX.sym("[a,b,c]")

#! Create a row matrix 
print SX.sym("[a,b,c]",1,3)

#! Space can be used as a separator, too:
print SX.sym("[a b c]")

#! Other brackets work fine as well:
print SX.sym("{a b c}")
print SX.sym("(a b c)")

#! Or you can omit the brackets entirely:
print SX.sym("a b c")

#! Create a matrix 
print SX.sym("[a b c;d e f]",2,3) 

#! The pythonic way to create a bunch of SX'es
#! ===========================================
#!
#! Thanks to tuple unpacking, you can write:
a,b,c = SX.sym("[a,b,c]")

#! This assigns three variables a, b and c:
print type(a), a
print type(b), b
print type(c), c

