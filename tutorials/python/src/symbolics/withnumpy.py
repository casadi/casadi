#! CasADi tutorial 1
#! ==================
#! This tutorial file explains the use of CasADi's SX in a python context.
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *

X=SXMatrix('X',2,2)
Y=SXMatrix('Y',2,2)
f=SXFunction ([X], [X])
print f.eval(Y)

a=SX("a")
b=SX("b")
c=SX("c")
d=SX("d")

A=array([[a,b],[c,d]])
B= SXMatrix(A)


