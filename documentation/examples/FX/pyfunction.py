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
def fac(f,nfwd,nadj,userdata):
  x = f.input()[0,0]
  y = 1
  for i in range(x):
    y*=(i+1)
  f.output().set(y)


#! Construct a PyFunction from it
c = PyFunction( fac, [sp_dense(1,1)], [sp_dense(1,1)] )
c.init()


#! Set the input
c.input().set(4)

#! Evaluate
c.evaluate(0,0)

print "4! = ", c.output().toScalar()

#! Using the function in a graph
#!==============================

x = MX("x")
[y] = c.call([x])

f = MXFunction([x],[y])
f.init()


c.input().set(5)
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

#! Jacobians
#!==============================
x = MX("x",3,1)
[y] = c.call([x])

f = MXFunction([x],[y])
f.init()


J = f.jacobian()
J.init()

c.setUserData({"foo": "bar"})

#! Our function will be called with forward seeds.
#! Of course, since we don't alter f.fwdSens, the result is meaningless
J.evaluate()




