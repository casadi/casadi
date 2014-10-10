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
#! CasADi tutorial 1
#! ==================
#! This tutorial file explains the use of CasADi's SXElement in a python context.
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *
#! Contructors & printing
#! --------------------------------------
#! Always provide a name to your symbols.
#! This name is not used for identification, only for printing.
a = SX.sym("z")
print type(a)
print a
#! You can explicitely create constant SX objects as follows
c = SX(5)
print c
#! Scalar algebra
#! --------------------------
#! Any operations on SX objects return SX objects.
x = SX.sym("x")
y = SX.sym("y")
c = x+y
print type(c)
print c
#! While you construct ever complex expressions, 
#! a graph of SX objects is created.
d = c*2 + x
print d
#! Note that, by itself, CasADi does very little expansions or simplifications of expressions.
#! Only simplifications that do not to introduce new nodes to the graph are allowed.
print d-x
print simplify(d-x)
print SX(5) + SX(7)
print 0*x + 0*y
print 1*x
#! SXElement objects are immutable entities.
#! The assignment and update operators are really creating new object instead of modifying the old ones.
print "object address before: %d" % id(d)
d = d - x
print d
print "object address after: %d" % id(d)
#! Consequently, updates and assignements don't have side effects for other SXElement objects
f = x + y
x *= 2
print x
print f
print x+y
#! Expression substitution
#! ------------------------------------
x=SX.sym("x")

y=x*x
print y
print substitute(y,x,SX.sym("w"))
print y
#! More operators
#! ------------------------------
#! Some familiar mathematical operations are supported that mimic the standard numpy functions:
#! sqrt sin cos tan arctan arcsin arccos exp log pow fabs floor ceil erf fmin fmax.
#! Using these functions require numpy to be imported.
y = sin(x**x)

x=SX.sym("x")
print type(x>0)
t = if_else(x>0,-10,10)
print t
#! Note that 'a<b' is treated as '!(a>=b)'


#! Conclusion
#! -------------------
#! We have seen how SXElement objects behave like symbolic objects.
#! They can be used to contruct expression trees.
#! 
#! To see how we can efficiently evaluate and differentiate these objects, jump on to the sxfunction tutorial...
