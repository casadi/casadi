from casadi import *

x = ssym("x")
f = sqrt(x)

F = SXFunction([x],[f])
F.setOption("verbose",True)
F.setOption("numeric_jacobian",False)
F.init()

J = F.jacobian()
J.setOption("verbose",True)
J.init()

J.setInput(0.738)
J.evaluate()
print "J.output(0) = ", J.output(0)
print "J.output(1) = ", J.output(1)

