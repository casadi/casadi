#! Symbolic substitution
#!======================
from casadi import *

#! Let's build a trivial symbolic SX graph
x = SX("x")
y = SX("y")
z_= x*y
z = z_+x 
print type(z), z

#! We need SXFuncion to manipulate the SX graph
f = SXFunction([[x,y]],[z])

#! We can substitute a leaf in the graph
w = SX("w")
q = f.eval([[w,y]])[0]
#! f.eval() returns a tuple with all outputs, we selected the first
print type(q), q
#! Note how q is now an SXMatrix

#! We can take a shortcut via substitute:
q = substitute([z],[x],[w])
print type(q), q

#! Note that substitution of non-symbolic SX nodes is not permitted:
#  substitute([z],[z_],[w])  This would throw an error
  
#! This is actually a restriction of SXFunction:
#  SXFunction([[z_,y]],[z])  This would throw an error
