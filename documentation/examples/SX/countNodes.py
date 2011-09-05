#! countNodes
#!======================
from casadi import *
from casadi.tools import *

#! Let's build a trivial symbolic SX graph
x = SX("x")
y = SX("y")
z = x*y+2*y
print countNodes(z), " nodes in ", z
dotdraw(z)

z += 4*z
print countNodes(z), " nodes in ", z
dotdraw(z)

z *= z+1
print countNodes(z), " nodes in ", z
dotdraw(z)


