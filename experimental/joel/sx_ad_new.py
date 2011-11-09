from casadi import *
from numpy import *

x = ssym("x",5)
f = sin(x) + x
print "f = ", f

fcn = SXFunction([x],[f])
fcn.init()

[f_recon,J] = fcn.jac([(0,-1),(0,0)])
print "f_recon = ", f_recon
print "J = ", J

x = SX("x")
y = SX("y")
z = SX("z")
w = SX("w")
f = [x, (x+(2*(y*y))), ((x+(2*(y*(y*y))))+(3*((z*z)*(z*z)))), w]
print "f = ", f

fcn = SXFunction([[x,y,z,w]],[f])
fcn.init()

[f_recon,J] = fcn.jac([(0,-1),(0,0)])
print "f_recon = ", f_recon
print "J = ", J



##q = SXMatrix(1,1,SX(5)) # segfault
#q = SXVector(1,SX(5)) # segfault
##q = DVector(1,5) # works
#d = q[:][0]
#print d
