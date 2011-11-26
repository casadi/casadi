from casadi import *
from casadi.tools import *

#! An SX graph
a = SX("a")
b = SX("b")

c = sin(a**5 + b)

c = c - b/ sqrt(fabs(c))
print c

dotdraw(c)

#! An MX graph
x = MX("x",sp_tril(2))
y = MX("y",sp_tril(2))

z = msym("z",4,2)

zz = x+y

f = MXFunction([z,y],[z+x[0],x-y])
f.setOption("name","magic")
f.init()

[z,z2] = f.call([vertcat([x,y]),zz.T])

z = z[:2,:] +x + cos(x) - sin(x) / tan(z2)

dotdraw(z)
