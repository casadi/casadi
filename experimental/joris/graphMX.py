from casadi import *
from casadi.tools import *

import __builtin__

x = msym("x",2,1)
y = msym("y")


z = (x+y)**4+sin(x)

dotsave(z,filename='MX1.pdf')


x = MX("x",sp_tril(2))

y = msym("y",2,2)

z = x+y

z = mul(vertcat([x,y]),z)


dotsave(z,filename='MX2.pdf')


x = MX("x",sp_tril(2))

y = msym("y",2,2)

z = x+y

f = MXFunction([x,y],[x+y,x-y])
f.setOption("name","my function")
f.init()

[z,z2] = f.call([x+y,x.T])

z = z+x + cos(x) - sin(x) / tan(z2)


dotsave(z,filename='MX3.pdf')


x = MX("x",sp_tril(2))

y = msym("y")

z = if_else(MX("y"),x**2,x**3)

dotsave(z,filename='MX4.pdf')

x = MX("x",4,1)

#z = norm_2(x)

dotsave(z,filename='MX5.pdf')


x = MX("x",3,3)
y = MX("y")

z = vec(x)

zz = z+y

zz = vertcat([zz,y])

d = [MX("d") for i in range(9)]
f = MXFunction(d + [MX("d",2,1)],[__builtin__.sum(i for i in d)])
f.init()

g = x[1:,1:]

[z] = f.call([z[0],z[1],z[2],z[3],x[4],x[5],zz[6],zz[7],zz[8],g[0:2]])

z = z +zz
dotsave(z,filename='MX5.pdf')
