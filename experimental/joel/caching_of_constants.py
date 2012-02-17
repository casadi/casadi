from casadi import *
a = ssym("a")
b = ((a+10)+10)+10
f = SXFunction([a],[b])
f.init()

print f

