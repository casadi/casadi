from casadi import *
from numpy import *

x=SX("x")
y=SX("y")
z=SX("z")
w=SX("w")

inp=[x,y,z,w]

out=SXMatrix(6,1)
out[0,0]=x
out[2,0]=x+2*y**2
out[4,0]=x+2*y**3+3*z**4
out[5,0]=w

reference = lambda x,y,z,w: array([[1,0,0,0],[0,0,0,0],[1,4*y,0,0],[0,0,0,0],[1,6*y**2,12*z**3,0],[0,0,0,1]])

f=SXFunction([inp],[out])
f.init()

print "inp=", inp
print "out=", out

Jf=f.jacobian(0,0)
print "We never get here"

Jf.init()

n=[1.2,2.3,7,4.6]
Jf.input().set(n)

Jf.evaluate()
J = reference(*n)

print "What it should be= ", J
print "What it is= ", Jf.output().toArray()
