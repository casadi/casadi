from casadi import *

#! SX
#! ===
#! Build up a huge expression
x = SX("x")
for i in range(5):
  x = sin(x)*x
  
#! With default printlimit
print x

SX.setMaxNumCallsInPrint(3)

print x


#! MX
#! ===
#! Build up a huge expression
x = MX("x")
for i in range(5):
  x = sin(x)*x
  
#! With default printlimit
print x

MX.setMaxNumCallsInPrint(3)

print x
