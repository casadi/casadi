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
#! IDAS integrator
#! =====================
#!
#! We solve a system
#!   $\dot{x}(t)=f(x(t),y(t),t)$ \n
#!   $0=g(x(t),y(t),t)$ \n
from casadi import *
from numpy import *
from pylab import *

#! We solve the following simple dae system that describes
#! the dynamics of a pendulum:
#! x' = u, y' = v, u' = lambda * x, v' =lambda * y - g
#!   s.t. x^2+y^2 = L
#!
#! We retain g and L as parameters
#! http://en.wikipedia.org/wiki/Differential_algebraic_equation#Examples
L = SX.sym("L")
g = SX.sym("g")

#! Time
t=SX.sym("t")

#! differential states
x=SX.sym("x")
y=SX.sym("y")
u=SX.sym("u")
v=SX.sym("v")

#! algebraic states
lambd=SX.sym("lambda")

#! All states and parameters
x_all = vertcat([x,u,y,v])
z_all = lambd
p_all = vertcat([L,g])

#! the initial state of the pendulum
P_ = [5,10] # parameters
X_ = [3,-1.0/3,4,1.0/4] # differential states
XDOT_ = [-1.0/3,1147.0/240,1.0/4,-653.0/180] # state derivatives
Z_ = [1147.0/720] # algebraic state

#! We construct the DAE system
ode = vertcat([u,lambd*x,v,lambd*y+g])
alg = x**2+y**2-L**2
f = SXFunction(daeIn(x=x_all,z=z_all,p=p_all,t=t),daeOut(ode=ode,alg=alg))
f.init()

#! Let's check we have consistent initial conditions:
f.setInput(P_,"p")
f.setInput(X_,"x")
f.setInput(Z_,"z")
f.evaluate()
print f.getOutput("ode") # This should be same as XDOT_
print f.getOutput("alg") # This should be all zeros


#! Let's check our jacobian $\frac{dg}{dy}$:
j = jacobian(alg,lambd)
print j
#! Note that the jacobian is not invertible: it is not of DAE-index 1
#!
#! This system is not solvable with idas, because it is of DAE-index 3.
#! It is impossible to lambda from the last element of the residual.

#! We create a DAE system solver
I = Integrator("idas", f)
I.setOption("calc_ic",False)
I.setOption("init_xdot",XDOT_)
  
I.init()

I.setInput(P_,"p")
I.setInput(X_,"x0")
I.setInput(Z_,"z0")

#! This system is not solvable with idas, because it is of DAE-index 3.
#! It is impossible obtain lambda from the last element of the residual.

try:
  I.evaluate()
except Exception as e:
  print e
  
#! We construct a reworked version od the DAE (index reduced), now it is DAE-index 1
ode = vertcat([u,lambd*x])
alg = vertcat([x**2+y**2-L**2, u*x+v*y,u**2-g*y+v**2+L**2*lambd])
x_all = vertcat([x,u])
z_all = vertcat([y,v,lambd])
f = SXFunction(daeIn(x=x_all,z=z_all,p=p_all,t=t),daeOut(ode=ode,alg=alg))
f.init()

#! the initial state of the pendulum
P_ = [5,10] # parameters
X_ = [3,-1.0/3] # differential states
XDOT_ = [-1.0/3,1147.0/240] # state derivatives
Z_ = [4,1.0/4,1147.0/720] # algebraic state

#! Let's check we have consistent initial conditions:
f.setInput(P_,"p")
f.setInput(X_,"x")
f.setInput(Z_,"z")
f.evaluate()
print f.getOutput("ode") # This should be the same as XDOT_
print f.getOutput("alg") # This should be all zeros

#! Let's check our jacobian:
j = f.jacobian("z","alg")
j.init()

j.setInput(P_,"p")
j.setInput(X_,"x")
j.setInput(Z_,"z")
j.evaluate()
print array(j.getOutput())
#! $\frac{dg}{dy}$ is invertible this time.

#! We create a DAE system solver
I = Integrator("idas", f)

I.setOption('t0',0)
I.setOption('tf',1)
I.setOption("init_xdot",XDOT_)
I.init()
  
I.setInput(P_,"p")
I.setInput(X_,"x0")
I.setInput(Z_,"z0")
I.evaluate()

print I.getOutput("xf")


#! Possible problems
#! ==================

#! If you would initialize with:
P_ = [5,10] # parameters
X_ = [5,0]  # states

#! You will get an error:
I.setInput(P_,"p")
I.setInput(X_,"x0")
I.setInput(Z_,"z0")
try:
  I.evaluate()
except Exception as e:
  print e 

#! Although this initialisation is consistent,
#! it coincides with a singular point.


 
