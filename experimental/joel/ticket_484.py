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
print "J.getOutput(0) = ", J.getOutput(0)
print "J.getOutput(1) = ", J.getOutput(1)

