from casadi import *
from casadi.tools import *

x = msym("x",2,1)
y = msym("y")


z = (x+y)**4+sin(x)

dotsave(z,filename='MX1.pdf')

z = mul(x.T,x)

dotsave(z,filename='MX2.pdf')
