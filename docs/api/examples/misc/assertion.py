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

f = MXFunction([x],[z])
f.init()

f.setInput(5)
f.evaluate()

print f.output()

#! For some mysterious reason we get NaN here

#! Next, we add an assertion: 
y = y.attachAssert(y>0,"bummer") # Add assertion here
z = sqrt(y)

f = MXFunction([x],[z])
f.init()

f.setInput(5)

try:
  f.evaluate()
except Exception as e:
  print "An exception was raised here:"
  print e


#! You can combine this with CustomFunction to do powerful assertions
@pyevaluate
def dummy(f):
  import numpy
  x = f.getInput()
  m = max(numpy.real(numpy.linalg.eig(blockcat([[x,-1],[-1,2]]))[0]))
  print "m=",m
  f.setOutput(int(m>2))


foo = CustomFunction(dummy, [x.sparsity()], [Sparsity.dense(1,1)] )
foo.init()

y = sin(x)

y = y.attachAssert(foo.call([y])[0],"you are in trouble") # Add assertion here
z = sqrt(y)

f = MXFunction([x],[z])
f.init()

f.setInput(5)

f.evaluate()




