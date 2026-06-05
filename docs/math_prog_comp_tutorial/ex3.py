# Python
import casadi as ca

# Formulate the ODE
x=ca.SX.sym('x',2)
p=ca.SX.sym('p')
z=1-x[1]**2
f=ca.vertcat(z*x[0]-x[1]+p,\
          x[0])
dae=dict(x=x,p=p,ode=f)

# Create solver instance
op=dict(t0=0,tf=1)
F=ca.integrator('F',\
      'cvodes',dae,op)

# Solve the problem
r=F(x0=[0,1],p=0.1)
print(r['xf'])

# Create Jacobian function
D=F.factory('D',\
 ['x0','p'],['jac:xf:x0'])

# Solve the problem
r=D(x0=[0,1],p=0.1)
print(r['jac_xf_x0'])
