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

#! CasADi provides a mechanism to add assertions in an MX expression graph
#! This can be useful to debug yor code, e.g. debugging why the end-result of a computation yields NaN

#! Consider this example:
x = MX.sym("x")
y = sin(x)
z = sqrt(y)

f = Function("f", [x], [z])

z0 = f(5)

print(z0)

#! For some mysterious reason we get NaN here

#! Next, we add an assertion:
y = y.attachAssert(y>0, "bummer") # Add assertion here
z = sqrt(y)

f = Function("f", [x],[z])

try:
  z0 = f(5)
except Exception as e:
  print("An exception was raised here:")
  print(e)


#! You can combine this with Callback to do powerful assertions
class Dummy(Callback):
  def __init__(self, name, opts={}):
    Callback.__init__(self)
    self.construct(name, opts)
  def get_n_in(self): return 1
  def get_n_out(self): return 1
  def eval(self, arg):
    import numpy
    x = arg[0]
    m = max(numpy.real(numpy.linalg.eig(blockcat([[x,-1],[-1,2]]))[0]))
    print("m=",m)
    return [int(m>2)]

foo = Dummy("foo")

y = sin(x)

y = y.attachAssert(foo(y), "you are in trouble") # Add assertion here
z = sqrt(y)

f = Function("f", [x],[z])

z0 = f(5)
