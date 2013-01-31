from casadi import *
a = msym("A",3,4)
b = msym("b",4,1)
c = mul(a,b)
c = 2 + c
f = MXFunction([a,b],[c])
f.init()
f.generateCode("f_mx.c")
