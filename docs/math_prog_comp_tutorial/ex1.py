# Python
from casadi import *

# Symbolic representation
x=SX.sym('x')
y=SX.sym('y')
z=y-(1-x)**2
f=x**2+100*z**2
P=dict(x=vertcat(x,y),f=f)

# Create solver instance
F=nlpsol('F','ipopt',P)

# Solve the problem
r=F(x0=[2.5,3.0])
print(r['x'])
