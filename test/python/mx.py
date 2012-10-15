from casadi import *

n = 3
m = 4
import numpy
numpy.random.seed(42)
A_ = numpy.random.random((m,n))
A = MX("A",m,n)
C_ = numpy.random.random((m,m))
C = MX("C",m,m)
D_ = numpy.random.random((m,n))
D = MX("D",m,n)

x_ = numpy.random.random((n,1))
x = MX("x",n,1)

Axb = mul(A,x)
Dxe = mul(D,x)
a = mul(mul(trans(Axb),C),Dxe)

f = MXFunction([x,A,C,D],[a])

for mode in ["forward", "reverse"]:
  f.init()
  J = f.jacobian()
  J.init()
  J.input(0).set(x_)
  J.input(1).set(A_)
  J.input(2).set(C_)
  J.input(3).set(D_)
  J.evaluate()
  
