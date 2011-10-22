from casadi import * 
from numpy import * 

x = SX("x")
y = SX("y")
z = SX("z")

var = SXMatrix([x,y,z])
dd = SXMatrix([5,sin(x),x+y])

print var
print dd

dd_rep = substituteInPlace(var,dd,True)
print dd_rep
