from casadi import *
x = ssym("x",3,4)
f1 = SXFunction([x],[sin(x)])
f1.init()

a = msym("A",3,4)
b = msym("b",4,1)
[z] = f1.call([a])
c = mul(z,b)
c = 2 + c
f = MXFunction([a,b],[c])
f.init()
f.generateCode("f_mx.c")
