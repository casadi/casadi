#! IDAS integrator
#! =====================
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
L = SX("L")
g = SX("g")

#! Time
t=SX("t")

#! differential states
x=SX("x")
y=SX("y")
u=SX("u")
v=SX("v")

dx=SX("dx")
dy=SX("dy")
du=SX("du")
dv=SX("dv")

#! algebraic states
lambd=SX("lambda")
dlambd=SX("dlambda")

#! the initial state of the pendulum
P_ = [5,10] # parameters
Y_ = [3,-1.0/3,4,1.0/4,1147.0/720] # states
YDOT_ = [-1.0/3,1147.0/240,1.0/4,-653.0/180,0] # state derivatives

#! We construct the DAE system
res = [dx-u,du-lambd*x,dy-v,dv-lambd*y+g,x**2+y**2-L**2]
f = SXFunction({DAE_T: t, DAE_P: [L,g], DAE_Y: [x,u,y,v,lambd],DAE_YDOT: [dx,du,dy,dv,dlambd]},[res])
f.init()

#! Let's check we have consistent initial conditions:
f.input(DAE_P).set(P_)
f.input(DAE_Y).set(Y_)
f.input(DAE_YDOT).set(YDOT_)
f.evaluate()
print f.output() # This should be all zeros


#! Let's check our jacobian:
j = Jacobian(f,DAE_Y,0)
j.init()

j.input(DAE_P).set(P_)
j.input(DAE_Y).set(Y_)
j.input(DAE_YDOT).set(YDOT_)
j.evaluate()
print array(j.output())

#! We create a DAE system solver
I = IdasIntegrator(f)
I.setOption("calc_ic",False)
I.setOption('is_differential',[1]*4 + [0])
  
I.init()

I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(Y_)
I.input(INTEGRATOR_XP0).set(YDOT_)
#! This system is not solvable with idas, because it is of DAE-index 3.
#! It is impossible to lambda from the last element of the residual.

try:
  I.evaluate()
except Exception as e:
  print e
  
#! We construct a reworked version od the DAE (index reduced), now it is DAE-index 1
res = [dx-u,du-lambd*x,x**2+y**2-L**2, u*x+v*y,u**2-g*y+v**2+L**2*lambd]
f = SXFunction({DAE_T: t, DAE_P: [L,g], DAE_Y: [x,u,y,v,lambd],DAE_YDOT: [dx,du,dy,dv,dlambd]},[res])
f.init()

#! Let's check we have consistent initial conditions:
f.input(DAE_P).set(P_)
f.input(DAE_Y).set(Y_)
f.input(DAE_YDOT).set(YDOT_)
f.evaluate()
print f.output() # This should be all zeros

#! Let's check our jacobian:
j = Jacobian(f,DAE_Y,0)
j.init()

j.input(DAE_P).set(P_)
j.input(DAE_Y).set(Y_)
j.input(DAE_YDOT).set(YDOT_)
j.evaluate()
print array(j.output())
#! The last three columns of the jacobian span a three-space, which allows us to solve. 

#! We create a DAE system solver
I = IdasIntegrator(f)
I.setOption('is_differential',[1]*2 + [0]*3)

I.setOption('t0',0)
I.setOption('tf',1)
I.init()
  
I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(Y_)
I.input(INTEGRATOR_XP0).set(YDOT_)
I.evaluate()

print I.output(INTEGRATOR_XF)


#! Possible problems
#! ==================

#! If you would initialize with:
P_ = [5,10] # parameters
Y_ = [5,0,0,0,0] # states
YDOT_ = [0,0,0,0,0] # state derivatives

#! You will get an error:
I.input(INTEGRATOR_P).set(P_)
I.input(INTEGRATOR_X0).set(Y_)
I.input(INTEGRATOR_XP0).set(YDOT_)
try:
  I.evaluate()
except Exception as e:
  print e 

#! Although this initialisation is consistent,
#! it coincides with a singular point.


 
