# Python
import casadi as ca

# Formulate the DAE
x=ca.SX.sym('x',2)
z=ca.SX.sym('z')
u=ca.SX.sym('u')
f=ca.vertcat(z*x[0]-x[1]+u,\
          x[0])
g=x[1]**2+z-1
h=x[0]**2+x[1]**2+u**2
dae=dict(x=x,p=u,ode=f,\
   z=z,alg=g,quad=h)

# Create solver instance
T = 10. # end time
N = 20 # discretization
op=dict(t0=0,tf=T/N)
F=ca.integrator('F',\
      'idas',dae,op)

# Empty NLP
w=[]; lbw=[]; ubw=[]
G=[]; J=0

# Initial conditions
Xk=ca.MX.sym('X0',2)
w+=[Xk]
lbw+=[0,1]
ubw+=[0,1]

for k in range(1,N+1):
  # Local control
  name='U'+str(k-1)
  Uk=ca.MX.sym(name)
  w+=[Uk]
  lbw+=[-1]
  ubw+=[ 1]

  # Call integrator
  Fk=F(x0=Xk,p=Uk)
  J+=Fk['qf']

  # Local state
  name='X'+str(k)
  Xk=ca.MX.sym(name,2)
  w+=[Xk]
  lbw+=[-.25,-ca.inf]
  ubw+=[ ca.inf, ca.inf]
  G+=[Fk['xf']-Xk]

# Create NLP solver
nlp=dict(f=J,\
     g=ca.vertcat(*G),\
     x=ca.vertcat(*w))
S=ca.nlpsol('S',\
     'blocksqp',nlp)

# Solve NLP
r=S(lbx=lbw,ubx=ubw,\
    x0=0,lbg=0,ubg=0)
print(r['x'])
