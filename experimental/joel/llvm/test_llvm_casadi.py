from casadi import *

# Construct a simple function
x = ssym("x")
f = sin(2*x)
F = SXFunction([x],[f])
F.setOption("just_in_time",True)
F.init()

