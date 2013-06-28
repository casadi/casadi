from casadi import *
from numpy import *
from os import system
import time
import sys

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
system("gcc -fPIC -shared f_mx.c  -o f_mx.so")

ef = ExternalFunction("./f_mx.so")
ef.init()

a_val = array([[1,2,3,3],[2,3,4,5],[3,4,5,6]])
b_val = array([7,6,5,4])

f.setInput(a_val,0);
f.setInput(b_val,1);
f.evaluate()

print f.getOutput()

ef.setInput(a_val,0);
ef.setInput(b_val,1);
ef.evaluate()

print f.getOutput()

