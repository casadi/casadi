# Python
import casadi as ca

# Formulate the NLP
x=ca.SX.sym('x')
y=ca.SX.sym('y')
z=ca.SX.sym('z')
f=x**2+100*z**2
g=z+(1-x)**2-y
P=dict(x=ca.vertcat(x,y,z),\
       f=f,g=g)

# Create solver instance
F=ca.nlpsol('F','ipopt',P)

# Solve the problem
r=F(x0=[2.5,3.0,0.75],\
    ubg=0, lbg=0)
print(r['x'])
