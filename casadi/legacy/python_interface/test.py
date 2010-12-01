# -*- coding: utf-8 -*-
from casadi_python import *
from numpy import *

x = MX('x')
y = MX('y')
f = x+sin(x*3)
print f

fcn = MXFunction(x,f)
print fcn


x = SXMatrix('x',2)
y = SXMatrix('y')
f = x+sin(y*3)
print f

fcn = SXFunction([x,y],f)
print fcn

# ode integrator
t = SXMatrix('t')
x = SXMatrix('x',2)
p = SXMatrix('p')

# ode rhs
xdot = x
ffcn = SXFunction([t,x,p],xdot)
print ffcn

# cvodes integrator
integrator = CVodesIntegrator(ffcn)
print integrator

