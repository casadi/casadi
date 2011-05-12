#! Code generation
#!======================
from casadi import *

#! Let's build a trivial symbolic SX graph
x = SX("x")
y = SX("y")
z = x*y+2*y
z += 4*z

#! An SXFunction is needed to inspect the graph
f = SXFunction([x,y],[z])

#! The default representation is just the name of the function
print f.__repr__()

#! A print statement will call __str__()
#! The result will look like a node-by-node tree evaluation
print f

#! The generateCode method will insert this node-by-node evaluation in exported C code
f.generateCode("generateCode.txt")

