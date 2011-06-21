#! expand
#!======================
from casadi import *
from numpy import *
import casadi as c

#! We construct a simple MX expression
x = MX("x",1,2)
y = MX("y",2,1)

z = c.prod(x,y)

#! Let's construct an MXfunction
f = MXFunction([x,y],[z])
f.init()

#! An MX graph is lazy in evaluation
print "Expression = ", f.outputMX()

#! We expand the MXFunction into an SXFunction
fSX = f.expand()

print "Expanded expression = ", fSX.outputSX()


#! Limitations
#! =============
#! Not all MX graphs can be expanded.
#! Here is an example of a situation where it will not work.
#!
#! Jacobian will perform numerical AD. We cannot expand an internal numerical algorithm.

j = Jacobian(f,0,0)
j.init()

J=MXFunction([x,y],j.call([x,y]))
J.init()

#! This function cannot be expanded.
try:
  J.expand()
except Exception as e:
  print e
  
#! Note that we could get the Jacobian of fSX as a workaround here.
