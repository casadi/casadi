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

#! A simple case of PyFunction
#!============================


#! We start with a python function with CFunctionWrapper-like arguments
#! The function calculates the factorial of its input
def fac(f,nfwd,nadj,userdata):
  x = f.input()[0,0]
  y = 1
  for i in range(x):
    y*=(i+1)
  f.setOutput(y)


#! Construct a PyFunction from it
c = PyFunction( fac, [sp_dense(1,1)], [sp_dense(1,1)] )
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

def dummy(f,nfwd,nadj,userdata):
  print "userdata: " , userdata
  if nfwd>0:
    print f.fwdSeed()
  
c = PyFunction(dummy, [sp_dense(3,1)], [sp_dense(3,1)] )
c.init()

c.setUserData(12345)
c.evaluate(0,0)

#! Of course, since we don't alter f.output, the result is meaningless
print c.output()


#! Providing sensitivities
#!==============================

#! This function calculates (x,y) -> (y+x**2,x*y)
#! Note that CasADi will in general request multiple seeding directions.
#! Our function must accomodate this.
def squares(f,nfwd,nadj,userdata):
  print "Called squares with :", (nfwd,nadj)
  x = f.input(0)[0]
  y = f.input(0)[1]

  f.setOutput(f.input(0)**2,0)
  f.setOutput(f.input(0)**2,0)
  
  for i in range(nfwd):
    xdot = f.fwdSeed(0,i)[0]
    ydot = f.fwdSeed(0,i)[1]
    f.setFwdSens([2*x*xdot+ydot,y*xdot+x*ydot],0,i)
    
  for i in range(nadj):
    xb = f.adjSeed(0,i)[0]
    yb = f.adjSeed(0,i)[1]
    f.setAdjSens([2*x*xb+y*yb,xb+x*yb],0,i)
    
c = PyFunction( squares, [sp_dense(2,1)], [sp_dense(2,1)] )
c.init()

#! Let's calculate the jacobian:

J = c.jacobian()
J.init()
J.setInput([3,5])
J.evaluate()

print J.output()

#! Forcing ad_mode is currently non-functional. See https://github.com/casadi/casadi/issues/614
c.setOption("ad_mode","reverse")
c.init()
J = c.jacobian()
J.init()
J.setInput([3,5])
J.evaluate()

print J.output()



