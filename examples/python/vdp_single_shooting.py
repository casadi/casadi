from casadi import *
from numpy import *
import matplotlib.pyplot as plt

nk = 20    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
t  = ssym("t")    # time
u  = ssym("u")    # control
x  = ssym("x",3)  # state
xp = ssym("xd",3) # state derivative

# ODE/DAE residual function
res = [(1 - x[1]*x[1])*x[0] - x[1] + u, \
       x[0], \
       x[0]*x[0] + x[1]*x[1] + u*u] - xp
f = SXFunction([t,x,u,xp],[res])

# Create an integrator
f_d = CVodesIntegrator(f)
f_d.setOption("abstol",1e-8) # tolerance
f_d.setOption("reltol",1e-8) # tolerance
f_d.setOption("steps_per_checkpoint",1000)
f_d.setOption("tf",tf/nk) # final time
f_d.init()

# All controls (use matrix graph)
U = msym("U",nk) # nk-by-1 symbolic variable

# The initial state (x_0=0, x_1=1, x_2=0)
X  = msym([0,1,0])

# State derivative (only relevant for DAEs)
Xp = msym([0,0,0])

# Build a graph of integrator calls
for k in range(nk):
  [X,Xp] = f_d.call([X,U[k],Xp])
  
# Objective function: x_2(T)
F = MXFunction([U],[X[2]])

# Terminal constraints: x_0(T)=x_1(T)=0
X_01 = X[0:2] # first two components of X
G = MXFunction([U],[X_01])

# Allocate an NLP solver
solver = IpoptSolver(F,G)
solver.init()

# Set bounds and initial guess
solver.setInput(-0.75*ones(nk), NLP_LBX)
solver.setInput(1.0*ones(nk), NLP_UBX)
solver.setInput(zeros(nk),NLP_X_INIT)
solver.setInput(zeros(2),NLP_LBG)
solver.setInput(zeros(2),NLP_UBG)

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = array(solver.output(NLP_X_OPT))

# Get values at the beginning of each finite element
tgrid = linspace(0,10,nk+1)
tgrid_u = linspace(0,10,nk)

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
