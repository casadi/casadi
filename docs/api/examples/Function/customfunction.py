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

#! A simple case of CustomFunction
#!================================


#! We start with a python function with evaluate-like arguments
#! The function calculates the factorial of its input
#! We annotate the function with the pyfunction decorator

@pyfunction([Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
def fac(inputs):
  x = inputs[0]
  y = 1
  for i in range(x):
    y*=(i+1)
  return [y]
  
#! Just as other CasADi function, we have to initialize
fac.init()

#! Set the input
fac.setInput(4)

#! Evaluate
fac.evaluate()

print "4! = ", fac.getOutput()

#! Using the function in a graph
#!==============================

x = MX.sym("x")
[y] = fac.call([x])

f = MXFunction([x],[y])
f.init()


f.setInput(5)
f.evaluate()

print "5! = ", f.getOutput()


#! Passing userdata
#!==============================
#! A clean way to pass around user data is to use a class instead of a function
#!
#! The class-approach is a bit more advanced in the following points:
#! - There is no decorator upfront, you have to make a call to PyFunction.
#! - The output of evaluation is not returned. Rather, you are given the outputs as mutable DMatrices that you should 'set'.
#!   This is so to make a very efficient allocation-free implementation possible.

class Fun:
  def __init__(self,userdata):
    self.userdata = userdata
    # sin(x+3*y)
        
  def evaluate(self,(x,y),(z,)):
    print "userdata: " , self.userdata
    
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    z.set(z2)
            
c = PyFunction(Fun(12345), [Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)] )
c.init()

c.setInput(1.2,0)
c.setInput(1.5,1)
c.evaluate()

print c.getOutput(), sin(1.2+3*1.5)

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
  
#! There are several ways to provide sensitivities.
#! The easiest way is using the class-approach and providing a 'fwd' and/or 'adj' method:

class Fun:
  # sin(x+3*y)
        
  def evaluate(self,(x,y),(z,)):
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    z.set(z2)

  def fwd(self,(x,y),(z,),seeds,sens):
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    z.set(z2)
    
    for ((dx,dy),(dz,)) in zip(seeds,sens):
      dz0 = 3*dy
      dz1 = dx+dz0
      dz2 = cos(z1)*dz1
      dz.set(dz2)
  
  def adj(self,(x,y),(z,),seeds,sens):
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    z.set(z2)
    
    for ((z_bar,),(x_bar,y_bar)) in zip(seeds,sens):
      bx = 0
      by = 0
      bz1 = 0
      bz0 = 0
      
      bz2 = z_bar
      bz1 += bz2*cos(z1)
      bx+= bz1;bz0+= bz1
      by+= 3*bz0
      x_bar.set(bx)
      y_bar.set(by)

c = PyFunction(Fun(), [Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)] )
c.init()

J = c.jacobian(1) # jacobian w.r.t second argument
J.init()
J.setInput(1.2,0)
J.setInput(1.5,1)
J.init()
J.evaluate()

print J.getOutput(), cos(1.2+3*1.5)*3

  
#! Another way to provide sensitivities is by providing the full Jacobian, i.e. the Jacobian of all inputs strung together with respect to all outputs strung together.
@pyfunction([Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
def fun((x,y)):
  # sin(x+3*y)

  z0 = 3*y
  z1 = x+z0
  z2 = sin(z1)

  return [z2]

fun.init()

@pyfunction([Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,2),Sparsity.dense(1,1)])
def funjac((x,y)):

  # First output, the full Jacobian
  J = horzcat((cos(x+3*y),3*cos(x+3*y)))

  # Remaining outputs, the original function outputs
  z0 = 3*y
  z1 = x+z0
  z2 = sin(z1)

  return [J,z2]

funjac.init()
fun.setFullJacobian(funjac)

J = fun.jacobian(1)
J.init()
J.setInput(1.2,0)
J.setInput(1.5,1)
J.evaluate()

print J.getOutput(0), cos(1.2+3*1.5)*3


#! The most advanced way is to provide a 'getDerivative' method in the class-approach that takes the number of fwd seeds and adjoint seeds and returns an Function.

class Fun:
  # sin(x+3*y)
  
  def evaluate(self,(x,y),(z,)):
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    z.set(z2)
    
  def getDerivative(self,f,nfwd,nadj):
    inputs = [f.input(i).sparsity() for i in range(f.getNumInputs())]
    outputs = [f.output(i).sparsity() for i in range(f.getNumOutputs())]
    
    sself = self

    class Der:
       def evaluate(self,xy_andseeds,z_andseeds):  sself.evaluateDer(xy_andseeds,z_andseeds,nfwd,nadj)

    FunDer = PyFunction(Der(),inputs+inputs*nfwd+outputs*nadj,outputs*(nfwd+1)+inputs*nadj)
    return FunDer
    
  def evaluateDer(self,inputs,outputs,nfwd,nadj):
    # sin(x+3*y)
    
    num_in  =  2
    num_out =  1
    
    x = inputs[0]
    y = inputs[1]
    
    z0 = 3*y
    z1 = x+z0
    z2 = sin(z1)
    outputs[0].set(z2)
    
    for i in range(nfwd):
      dx = inputs[num_in + i*num_in+0]
      dy = inputs[num_in + i*num_in+1]
      
      dz0 = 3*dy
      dz1 = dx+dz0
      dz2 = cos(z1)*dz1
      
      outputs[num_out + i].set(dz2)
    
    for i in range(nadj):
      # Backwards sweep
      bx = 0
      by = 0
      bz1 = 0
      bz0 = 0
      
      bz2 = inputs[num_in + nfwd*num_in+i*num_out+0]
      bz1 += bz2*cos(z1)
      bx+= bz1;bz0+= bz1
      by+= 3*bz0
      outputs[num_out + nfwd*num_out + num_in*i+0].set(bx)
      outputs[num_out + nfwd*num_out + num_in*i+1].set(by)
    

Fun = PyFunction(Fun(),[Sparsity.dense(1,1),Sparsity.dense(1,1)], [Sparsity.dense(1,1)])
Fun.init()
J = Fun.jacobian(1)
J.init()
J.setInput(1.2,0)
J.setInput(1.5,1)
J.evaluate()

print J.getOutput(0), cos(1.2+3*1.5)*3

#! Note that this last approach allows to make return 'getDerivative' another PyFunction which in turn implements its own 'getDerivative' in order to provide higher prder sensitivities.

