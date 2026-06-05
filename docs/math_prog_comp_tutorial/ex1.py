# Python
import casadi as ca

# Symbolic representation
x=ca.SX.sym('x')
y=ca.SX.sym('y')
z=y-(1-x)**2
f=x**2+100*z**2
P=dict(x=ca.vertcat(x,y),f=f)

# Create solver instance
F=ca.nlpsol('F','ipopt',P)

# Solve the problem
r=F(x0=[2.5,3.0])
print(r['x'])
