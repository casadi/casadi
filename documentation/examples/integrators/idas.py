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
L = ssym("L")
g = ssym("g")

#! Time
t=ssym("t")

#! differential states
x=ssym("x")
y=ssym("y")
u=ssym("u")
v=ssym("v")

dx=ssym("dx")
dy=ssym("dy")
du=ssym("du")
dv=ssym("dv")

#! algebraic states
lambd=ssym("lambda")

#! All states and parameters
x_all = vertcat([x,u,y,v])
dx_all = vertcat([dx,du,dy,dv])
z_all = lambd
p_all = vertcat([L,g])

#! the initial state of the pendulum
P_ = [5,10] # parameters
X_ = [3,-1.0/3,4,1.0/4] # differential states
XDOT_ = [-1.0/3,1147.0/240,1.0/4,-653.0/180] # state derivatives
Z_ = [1147.0/720] # algebraic state

#! We construct the DAE system
ode = vertcat([dx-u,du-lambd*x,dy-v,dv-lambd*y+g])
alg = x**2+y**2-L**2
f = SXFunction(daeIn(x=x_all,z=z_all,p=p_all,t=t,xdot=dx_all),daeOut(ode=ode,alg=alg))
f.init()

#! Let's check we have consistent initial conditions:
f.input(DAE_P).set(P_)
f.input(DAE_X).set(X_)
f.input(DAE_Z).set(Z_)
f.input(DAE_XDOT).set(XDOT_)
f.evaluate()
print f.output(DAE_ODE) # This should be all zeros
print f.output(DAE_ALG) # This should be all zeros


#! Let's check our jacobian $\frac{dg}{dy}$:
j = jacobian(alg,lambd)
print j
#! Note that the jacobian is not invertible: it is not of DAE-index 1
#!
#! This system is not solvable with idas, because it is of DAE-index 3.
#! It is impossible to lambda from the last element of the residual.

#! We create a DAE system solver
I = IdasIntegrator(f)
I.setOption("calc_ic",False)
I.setOption("init_z",Z_)
I.setOption("init_xdot",XDOT_)
  
I.init()

I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(X_)

#! This system is not solvable with idas, because it is of DAE-index 3.
#! It is impossible obtain lambda from the last element of the residual.

try:
  I.evaluate()
except Exception as e:
  print e
  
#! We construct a reworked version od the DAE (index reduced), now it is DAE-index 1
ode = vertcat([dx-u,du-lambd*x])
alg = vertcat([x**2+y**2-L**2, u*x+v*y,u**2-g*y+v**2+L**2*lambd])
x_all = vertcat([x,u])
dx_all = vertcat([dx,du])
z_all = vertcat([y,v,lambd])
f = SXFunction(daeIn(x=x_all,z=z_all,p=p_all,t=t,xdot=dx_all),daeOut(ode=ode,alg=alg))
f.init()

#! the initial state of the pendulum
P_ = [5,10] # parameters
X_ = [3,-1.0/3] # differential states
XDOT_ = [-1.0/3,1147.0/240] # state derivatives
Z_ = [4,1.0/4,1147.0/720] # algebraic state

#! Let's check we have consistent initial conditions:
f.input(DAE_P).set(P_)
f.input(DAE_X).set(X_)
f.input(DAE_Z).set(Z_)
f.input(DAE_XDOT).set(XDOT_)
f.evaluate()
print f.output(DAE_ODE) # This should be all zeros
print f.output(DAE_ALG) # This should be all zeros

#! Let's check our jacobian:
j = f.jacobian(DAE_Z,DAE_ALG)
j.init()

j.input(DAE_P).set(P_)
j.input(DAE_X).set(X_)
j.input(DAE_Z).set(Z_)
j.input(DAE_XDOT).set(XDOT_)
j.evaluate()
print array(j.output())
#! $\frac{dg}{dy}$ is invertible this time.

#! We create a DAE system solver
I = IdasIntegrator(f)

I.setOption('t0',0)
I.setOption('tf',1)
I.setOption("init_z",Z_)
I.setOption("init_xdot",XDOT_)
I.init()
  
I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(X_)
I.evaluate()

print "output = ", I.output(INTEGRATOR_XF)
print "z = ", I.output(INTEGRATOR_ZF)

#! Possible problems
#! ==================

#! If you would initialize with:
P_ = [5,10] # parameters
X_ = [5,0]  # states

#! You will get an error:
I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(X_)
try:
  I.evaluate()
except Exception as e:
  print e 

#! Although this initialisation is consistent,
#! it coincides with a singular point.


 
