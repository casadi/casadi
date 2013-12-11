#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#! A simple case of CustomFunction
#!================================


#! We start with a python function with evaluate-like arguments
#! The function calculates the factorial of its input
#! We annotate the function with the pyevaluate decorator

@pyevaluate
def fac(f):
  x = f.getInput()[0,0]
  y = 1
  for i in range(x):
    y*=(i+1)
  f.setOutput(y)


#! Construct a CustomFunction from it
c = CustomFunction( fac, [sp_dense(1,1)], [sp_dense(1,1)] )
c.init()


#! Set the input
c.setInput(4)

#! Evaluate
c.evaluate()

print "4! = ", c.output().toScalar()

#! Using the function in a graph
#!==============================

x = MX("x")
[y] = c.call([x])

f = MXFunction([x],[y])
f.init()


c.setInput(5)
c.evaluate()

print "5! = ", c.output().toScalar()


#! Passing userdata
#!==============================
#! A clean way to pass around user data is to use a class
#! and make it callable py implementing __call__

@pyevaluate
class Dummy:
  def __init__(self,user_data):
    self.userdata = user_data
  def __call__(self,f):
    print "userdata: " , self.userdata
  
c = CustomFunction(Dummy(12345), [sp_dense(3,1)], [sp_dense(3,1)] )
c.init()


c.evaluate()

#! Of course, since we don't alter f.output, the result is meaningless
print c.getOutput()

#! Providing sensitivities
#!==============================
#! 
#! By default, there is no way to calculate derivatives from a function given as a CustomFunction
#! Operations requiring derivative information therefore result in an error:

try:
  J = c.jacobian()
  J.init()
except Exception as e:
  print e
  
#! We can provide forward sensitivities by providing the full Jacobian, i.e. the Jacobian of all inputs with respect to all outputs.
@pyevaluate
def fun(f):
  # sin(x+3*y)
  
  x = f.input(0)
  y = f.input(1)
  
  z0 = 3*y
  z1 = x+z0
  z2 = sin(z1)

  f.setOutput(z2)

c = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)])
c.init()

def funjac(f):
  x = f.input(0)
  y = f.input(1)

  # First output, the full Jacobian
  J = horzcat((cos(x+3*y),3*cos(x+3*y)))
  f.setOutput(J,0)

  # Remaining outputs, the original function outputs
  z0 = 3*y
  z1 = x+z0
  z2 = sin(z1)

  f.setOutput(z2,1)

jc = CustomFunction(funjac, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,2),sp_dense(1,1)])
jc.init()
c.setFullJacobian(jc)

J = c.jacobian()
J.init()
J.setInput(0.3,0)
J.setInput(0.7,0)
J.evaluate()

print J.getOutput()
