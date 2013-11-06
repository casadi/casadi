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
def fac(f,nfwd,nadj):
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
c.evaluate(0,0)

print "4! = ", c.output().toScalar()

#! Using the function in a graph
#!==============================

x = MX("x")
[y] = c.call([x])

f = MXFunction([x],[y])
f.init()


c.setInput(5)
c.evaluate(0,0)

print "5! = ", c.output().toScalar()


#! Passing userdata
#!==============================
#! A clean way to pass around user data is to use a class
#! and make it callable py implementing __call__

@pyevaluate
class Dummy:
  def __init__(self,user_data):
    self.userdata = user_data
  def __call__(self,f,nfwd,nadj):
    print "userdata: " , self.userdata
    if nfwd>0:
      print f.getFwdSeed()
  
c = CustomFunction(Dummy(12345), [sp_dense(3,1)], [sp_dense(3,1)] )
c.init()


c.evaluate(0,0)

#! Of course, since we don't alter f.output, the result is meaningless
print c.getOutput()

#! Providing sensitivities
#!==============================
#! 
#! The arguments nfwd and nadj stand for number of forward derivatives and number of adjoint derivatives respectively.
#! By default, they will always be zero. Indeed, doing an action that would result in CasADi requesting non-zero nfwd, nadj results in an error:

try:
  J = c.jacobian()
  J.init()
except Exception as e:
  print e
  
#! We can provide forward senitivities as follows:
#! This function calculates (x,y) -> sin(x+3*y)

@pyevaluate
def fun(f,nfwd,nadj):
  # sin(x+3*y)
  print "Called fun with :", (nfwd,nadj)
  
  x = f.input(0)
  y = f.input(1)
  
  dx = f.fwdSeed(0)
  dy = f.fwdSeed(1)
  
  z0 = 3*y
  dz0 = 3*dy
  
  z1 = x+z0
  dz1 = dx+dz0
  
  z2 = sin(z1)
  dz2 = cos(z1)*dz1

  # Backwards sweep
  bx = 0
  by = 0
  bz1 = 0
  bz0 = 0
  
  bz2 = f.adjSeed(0)
  bz1 += bz2*cos(z1)
  bx+= bz1;bz0+= bz1
  by+= 3*bz0
  f.setAdjSens(bx,0)
  f.setAdjSens(by,1)
  
  f.setOutput(z2)
  f.setFwdSens(dz2)

c = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)] )
c.setOption("max_number_of_fwd_dir",1)
c.setOption("max_number_of_adj_dir",1)
c.init()

#! Now this function is capable of delivering numeric sensitivity info.
#! It can be used to calculate a gradient for example:

J = c.jacobian()
J.init()

J.setInput(0.3,0)
J.setInput(0.7,0)
J.evaluate()

print J.output()
      
#! Forcing reverse mode:
c.setOption("ad_mode","reverse")
c.init()
J = c.jacobian()
J.init()
J.setInput(0.3,0)
J.setInput(0.7,0)
J.evaluate()

print J.getOutput()

#! We can also provide a function that returns the jacobian, instead of specifying the seeds

@pyevaluate
def fun(f,nfwd,nadj):
  # sin(x+3*y)
  print "Called fun with :", (nfwd,nadj)
  
  x = f.input(0)
  y = f.input(1)
  
  z0 = 3*y
  z1 = x+z0
  z2 = sin(z1)

  f.setOutput(z2)

@jacobiangenerator
def funjac(f,iind,oind):
  # sin(x+3*y)
  print "Called jacobian with :", (iind,oind)
  
  x = ssym("x")
  y = ssym("y")
  
  if iind == 0:
    J = cos(x+3*y)
  elif iind==1:
    J = 3*cos(x+3*y)
    
  f = SXFunction([x,y],[J])
  f.init()
  
  return f
  
c = CustomFunction(fun, [sp_dense(1,1),sp_dense(1,1)], [sp_dense(1,1)] )
c.setOption("jacobian_generator",funjac)
c.init()

J = c.jacobian()
J.init()
J.setInput(0.3,0)
J.setInput(0.7,0)
J.evaluate()

print J.getOutput()
