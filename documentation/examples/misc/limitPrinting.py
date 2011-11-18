from casadi import *

#! SX
#! ===
#! Build up a huge expression
x = SX("x")
for i in range(5):
  x = sin(x)*x
  
#! With default printlimit
print x

cvar.SX_max_num_calls_in_print = 3

print x


#! MX
#! ===
#! Build up a huge expression
x = MX("x")
for i in range(5):
  x = sin(x)*x
  
#! With default printlimit
print x

cvar.MX_max_num_calls_in_print = 3

print x
