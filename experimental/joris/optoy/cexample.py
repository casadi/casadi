from casadi import *
from coptoy import *

# Example 1: SDP problem
x = var()
y = var()

minimize(2*x+y,[blockcat([[-x,2],[2,-y]])<=0])

print x.sol, y.sol # sol: sqrt(2) , 2*sqrt(2)

# Example 2: linear system
x0 = var(lb=0)
x1 = var(lb=0)

minimize(2*x0+x1*3,[x0+x1>=1])  

print x0.sol, x1.sol # 1 0


