from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

# Declare variables (use simple, efficient DAG)
t = SX("t") # time
x=SX("x"); y=SX("y"); u=SX("u"); L=SX("cost")

# ODE right hand side function
f = [(1 - y*y)*x - y + u, x, x*x + y*y + u*u]
rhs = SXFunction([[t],[x,y,L],[u]],[f])

# Create an integrator (CVodes)
I = CVodesIntegrator(rhs)
I.setOption("abstol",1e-8) # abs. tolerance
I.setOption("reltol",1e-8) # rel. tolerance
I.setOption("steps_per_checkpoint",50)
I.setOption("stop_at_end",True)
I.init()

# Final time (fixed)
tf = 20.0

# Control bounds
u_min = -0.75
u_max = 1.0
u_init = 0.0

# State bounds and initial guess
x_min = -NP.inf;  x_max =  NP.inf;  x_init = 0
y_min = -NP.inf;  y_max =  NP.inf;  y_init = 0
L_min = -NP.inf;  L_max =  NP.inf;  L_init = 0

# State bounds at the final time
xf_min = 0;        xf_max =  0
yf_min = 0;        yf_max =  0
Lf_min = -NP.inf;  Lf_max =  NP.inf

# Numboer of shooting nodes
NS = 50

# Number of discretized controls
NU = NS

# Number of discretized states
NX = 3*NS

# Declare variable vector
NV = NU+NX
V = MX("V",NV)

# Disretized control
U = []
U_min = []
U_max = []
U_init = []
for i in range(NS):
  U.append(V[i])
  U_min.append(u_min)
  U_max.append(u_max)
  U_init.append(u_init)

# Disretized state
X = []
X_min = []
X_max = []
X_init = []
for i in range(NS):
  X += [V[NU+i*3 : NU+(i+1)*3]]
  X_init += [x_init, y_init, L_init]
  if i==NS-1:
    # final time
    X_min += [xf_min, yf_min, Lf_min]
    X_max += [xf_max, yf_max, Lf_max]
  else:
    # interor time
    X_min += [x_min, y_min, L_min]
    X_max += [x_max, y_max, L_max]
  
# Beginning of each shooting interval (shift time horizon)
T0 = NS*[MX(0)]

# End of each shooting interval (shift time horizon)
TF = NS*[MX(tf/NU)]

# The initial state (x=0, y=1, L=0)
X0  = MX([0,1,0])

# State derivative (not used)
XP = MX()

# Algebraic state (not used)
Z = MX()

# Constraint function with upper and lower bounds
g = []
g_min = []
g_max = []

# Build up a graph of integrator calls
for k in range(NS):
  # call the integrator
  [XF,XP,Z] = I.call([T0[k],TF[k],X0,U[k],XP,Z])
  
  # append continuity constraints
  g.append(X[k] - XF)
  g_min.append(NP.zeros(XF.numel()))
  g_max.append(NP.zeros(XF.numel()))
  
  # Update initial state for next interval
  X0 = X[k]

# State at the final time
XF = X[NS-1]

# Objective function: L(T)
F = MXFunction([V],[XF[2]])

# Terminal constraints: 0<=[x(T);y(T)]<=0
G = MXFunction([V],[vertcat(g)])

## Jacobian of a block
#Jxx = I.jacobian(INTEGRATOR_X0,INTEGRATOR_XF)
#Jxp = I.jacobian(INTEGRATOR_P, INTEGRATOR_XF)

## Build up a graph of integrator calls
#X0  = MX([0,1,0])
#for k in range(NS):
  ## call the integrator
  #Jk = Jxx.call([T0[k],TF[k],X0,U[k],XP,Z])
  #X0 = X[k]

## Construct the Jacobian (quickfix)

solver = IpoptSolver(F,G)
solver.setOption("tol",1e-5)
solver.setOption("hessian_approximation", "limited-memory")
solver.setOption("max_iter",100)
solver.setOption("linear_solver","ma57")
#solver.setOption("verbose",True)
solver.init()

# Set bounds and initial guess
solver.setInput(U_min  + X_min,  NLP_LBX)
solver.setInput(U_max  + X_max,  NLP_UBX)
solver.setInput(U_init + X_init, NLP_X_INIT)
solver.setInput(NP.concatenate(g_min),NLP_LBG)
solver.setInput(NP.concatenate(g_max),NLP_UBG)

# Solve the problem
solver.solve()

# Retrieve the solution
v_opt = solver.output(NLP_X_OPT)
u_opt = v_opt[0:NU]
X_opt = v_opt[NU:]
x_opt = [X_opt[k] for k in range(0,3*NS,3)]
y_opt = [X_opt[k] for k in range(1,3*NS,3)]
l_opt = [X_opt[k] for k in range(2,3*NS,3)]

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(u_opt)
plt.title("u_opt")

plt.figure(2)
plt.clf()
plt.plot(x_opt)
plt.title("x_opt")

plt.figure(3)
plt.clf()
plt.plot(y_opt)
plt.title("y_opt")

plt.figure(4)
plt.clf()
plt.plot(l_opt)
plt.title("l_opt")


#plt.show()

