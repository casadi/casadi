from casadi import *

x = MX.sym("x",2)

y = vertsplit(x,[0,1,2])[1]

f = MXFunction([x],[y])
f.init()

H = f.hessian()
H.init()
