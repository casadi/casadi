from casadi import *
from casadi.tools import *

x = SX("x")
y = SX("y")

dotsave(x**8,filename='SX0.pdf')

z = jacobian(sqrt(x**2+y**2),x)
print "done"

dotsave(z,filename='SX1.pdf')

z = jacobian(sqrt(x**2+y**2),[x,y])

dotsave(z,filename='SX2.pdf')

z = SXMatrix(4,5)
z[:,1] = x
z[1,:] = x**2+y**2
z[3,3] = 66.6
z[3,2] = cos(y)

dotsave(z,filename='SX3.pdf')

