#! jacSparsity
#!======================
from casadi import *
from numpy import *
import casadi as c
from pylab import spy, show

#! We construct a simple SX expression
x = ssym("x",40)
y = x[:-2]-2*x[1:-1]+x[2:]

#! Let's see what the first 5 entries of y look like
print y[:5]

#! Next, we construct a function
f = SXFunction([x],[y])
f.init()

#! And we visualize the sparsity of the jacobian
spy(f.jacSparsity())

show()

