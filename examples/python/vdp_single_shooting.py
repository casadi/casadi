from casadi import *
from numpy import *
import matplotlib.pyplot as plt

# Number of shooting nodes
NS = 20

# Declare variables (use simple, efficient DAG)
t = SX("t") # time
x=SX("x"); y=SX("y"); u=SX("u"); L=SX("cost")

# ODE right hand side function
f = [(1 - y*y)*x - y + u, x, x*x + y*y + u*u]
rhs = SXFunction([[t],[x,y,L],[u],[]],[f])

# Create an integrator (CVodes)
I = CVodesIntegrator(rhs)
I.setOption("abstol",1e-8) # abs. tolerance
I.setOption("reltol",1e-8) # rel. tolerance
I.setOption("steps_per_checkpoint",1000)
I.setOption("fsens_err_con",False)
I.setOption("stop_at_end",False)
I.setOption("t0",0.0)
I.setOption("tf",10.0/NS)
I.init()

# All controls (use complex, general DAG)
U = MX("U",NS)

# The initial state (x=0, y=1, L=0)
X  = MX([0,1,0])

# State derivative (not used)
XP = MX()

# Build up a graph of integrator calls
for k in range(NS):
  # Call the integrator
  [X,XP] = I.call([X,U[k],XP])
  
# Objective function: L(T)
cost = X[2]
F = MXFunction([U],[cost])

# Terminal constraints: 0<=[x(T);y(T)]<=0
eq = X[0:2]
G = MXFunction([U],[eq])

## Lagrange multipliers
#lam = MX("lam",eq.size1())

## Objective scaling factor
#sigma = MX("sigma")

## Lagrange function
#ll = sigma*obj + inner_prod(lam,eq)
#L = MXFunction([U,lam,sigma],[ll])
#L.init()

## Gradient of the Lagrangian
#lg = trans(L.grad()[0])
#GL = MXFunction([U,lam,sigma],[lg])
#GL.init()

## Hessian of the Lagrangian
#H = GL.jacobian()

# Allocate NLP solver
solver = IpoptSolver(F,G)
solver.setOption("hessian_approximation", \
                "limited-memory")
  
# Initialize the NLP solver
solver.init()

# Set bounds and initial guess
solver.setInput(NS*[-0.75], NLP_LBX)
solver.setInput(NS*[1.0],NLP_UBX)
solver.setInput(NS*[0.0],NLP_X_INIT)
solver.setInput([0,0],NLP_LBG)
solver.setInput([0,0],NLP_UBG)

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = array(solver.output(NLP_X_OPT))

# Get values at the beginning of each finite element
tgrid = linspace(0,10,NS+1)
tgrid_u = linspace(0,10,NS)

# Plot the results
plt.figure(1)
plt.clf()
#plt.plot(tgrid,x_opt,'--')
#plt.plot(tgrid,y_opt,'-')
plt.plot(tgrid_u,u_opt,'-.')
plt.title("Van der Pol optimization - single shooting")
plt.xlabel('time')
#plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.legend(['u trajectory'])
plt.grid()
plt.show()
