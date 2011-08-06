# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt

# Declare variables (use simple, efficient DAG)
t = SX("t") # time
xx=SX("x"); yy=SX("y"); uu=SX("u"); ll=SX("cost")
x = [xx,yy,ll]
u = [uu]

# Right hand side of the ODE
xdot = [(1 - yy*yy)*xx - yy + uu, xx, xx*xx + yy*yy + uu*uu]
ffcn = SXFunction([[t],x,u],[xdot])

# Objective function (meyer term)
mfcn = SXFunction([[t],x,u],[[ll]])

# Nonlinear constraint function
cfcn = FX() # Not used in this example

# Final time (fixed)
tf = 20.0

# Control bounds
uu_min = -0.75
uu_max = 1.0
uu_init = 0.0

u_lb = NP.array([uu_min])
u_ub = NP.array([uu_max])
u_init = NP.array([uu_init])

# State bounds and initial guess
xx_min = -NP.inf;  xx_max =  NP.inf;  xx_init = 0
yy_min = -NP.inf;  yy_max =  NP.inf;  yy_init = 0
ll_min = -NP.inf;  ll_max =  NP.inf;  ll_init = 0

# State bounds at the final time
xxf_min = 0;        xxf_max =  0
yyf_min = 0;        yyf_max =  0
llf_min = -NP.inf;  llf_max =  NP.inf

# Bounds on the states
x_lb = NP.array([xx_min,yy_min,ll_min])
x_ub = NP.array([xx_max,yy_max,ll_max])
x_init = NP.array([xx_init,yy_init,ll_init])
xf_lb = NP.array([xxf_min,yyf_min,llf_min])
xf_ub = NP.array([xxf_max,yyf_max,llf_max])

# Initialize functions
ffcn.init()  
mfcn.init()
#cfcn.init()

# Dimensions
nx = len(x)
nu = len(u)
nc = 0

# Legendre collocation points
legendre_points1 = [0,0.500000]
legendre_points2 = [0,0.211325,0.788675]
legendre_points3 = [0,0.112702,0.500000,0.887298]
legendre_points4 = [0,0.069432,0.330009,0.669991,0.930568]
legendre_points5 = [0,0.046910,0.230765,0.500000,0.769235,0.953090]
legendre_points = [0,legendre_points1,legendre_points2,legendre_points3,legendre_points4,legendre_points5]

# Radau collocation points
radau_points1 = [0,1.000000]
radau_points2 = [0,0.333333,1.000000]
radau_points3 = [0,0.155051,0.644949,1.000000]
radau_points4 = [0,0.088588,0.409467,0.787659,1.000000]
radau_points5 = [0,0.057104,0.276843,0.583590,0.860240,1.000000]
radau_points = [0,radau_points1,radau_points2,radau_points3,radau_points4,radau_points5]

# Type of collocation points
LEGENDRE = 0
RADAU = 1
collocation_points = [legendre_points,radau_points]

# Degree of interpolating polynomial
K = 3

# Number of finite elements
N = 50

# Radau collocation points
cp = RADAU

# Size of the finite elements
h = tf/N

# Coefficients of the collocation equation
C = (K+1) * [[]]

# Coefficients of the continuity equation
D = (K+1) * [[]]

# Collocation point
tau = SX("tau")
  
# Roots
tau_root = collocation_points[cp][K]

for j in range(K+1):
  # Lagrange polynomials
  L = SX(1)
  for k in range(K+1):
    if k != j:
      L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k])
  
    # Lagrange polynomial function
    lfcn = SXFunction([[tau]],[[L]])
    lfcn.init()
    
    # Get the coefficients of the continuity equation
    lfcn.setInput(1.0)
    lfcn.evaluate()
    D[j] = lfcn.output()[0]

    # Get the coefficients of the collocation equation
    C[j] = (K+1) * [[]]
    for k in range(K+1):
      lfcn.setInput(tau_root[k])
      lfcn.setFwdSeed(1.0)
      lfcn.evaluate(1,0)
      C[j][k] = lfcn.fwdSens()[0]

# Collocated times
T = N * [[]]
for i in range(N):
  T[i] = (K+1) * [[]]
  for j in range(K+1):
    T[i][j] = h*(i + collocation_points[cp][K][j])

# Total number of variables
NX = N*(K+1)*nx
NU = N*nu
NXF = nx
NV = NX+NU+NXF

# NLP variable vector
V = MX("V",NV)

# Collocated states
X = []

# Parametrized control (piecewice constant)
U = []
  
# All variables with bounds and initial guess
vars_lb = NP.zeros(NV)
vars_ub = NP.zeros(NV)
vars_init = NP.zeros(NV)
vi = 0

# Loop over the finite elements
for i in range(N):  
  # Collocated states
  Xi = []
  for j in range(K+1):
    Xi.append(V[vi:vi+nx])
    vars_init[vi:vi+nx] = x_init
    if i==0 and j==0:
      vars_lb[vi:vi+nx] = NP.array([0,1,0])
      vars_ub[vi:vi+nx] = NP.array([0,1,0])
    else:
      vars_lb[vi:vi+nx] = x_lb
      vars_ub[vi:vi+nx] = x_ub
    vi += nx
        
  X.append(Xi)
  
  # Parametrized controls
  U.append(V[vi:vi+nu])
  vars_lb[vi:vi+nu] = u_lb
  vars_ub[vi:vi+nu] = u_ub
  vars_init[vi:vi+nu] = u_init
  vi += nu
  
# State at end time
XF = V[vi:vi+nx]
vars_lb[vi:vi+nx] = xf_lb
vars_ub[vi:vi+nx] = xf_ub
vars_init[vi:vi+nx] = x_init
vi += nx
  
# Constraint function for the NLP
g = []
lbg = []
ubg = []

for i in range(N):
  for k in range(1,K+1):
    # augmented state vector
    y_ik = [MX(T[i][k]), X[i][k], U[i]]
    
    # Add collocation equations to NLP
    [temp] = ffcn.call(y_ik)
    for j in range(temp.size()):
      temp[j] *= h
    for j in range (K+1):
      for l in range(temp.size()):
        temp[l] -= X[i][j][l]*C[j][k]
      
    g += [temp]
    lbg.append(NP.zeros(nx)) # equality constraints
    ubg.append(NP.zeros(nx)) # equality constraints
      
    # Add nonlinear constraints, if any
    if nc>0:
      g += cfcn.call(y_ik)
      lbg.append(c_lb)
      ubg.append(c_ub)

  # Add continuity equation to NLP
  if i<N-1:
    temp = MX(X[i+1][0])
  else:
    temp = MX(XF)
    
  for j in range(K+1):
    for l in range(nx):
      temp[l] -= D[j]*X[i][j][l]

  g += [temp]
  lbg.append(NP.zeros(nx))
  ubg.append(NP.zeros(nx))
  
# Variable vector (SX)
V_sx = symbolic("V",NV)
  
# Nonlinear constraint function
gg = vertcat(g)
gfcn_nlp_mx = MXFunction([V],[gg])

# Objective function of the NLP
y_f = [MX(T[N-1][K]),XF,U[N-1]]
[f] = mfcn.call(y_f)
ffcn_nlp_mx = MXFunction([V], [f])

# Expand constraint function as an SX graph (move to NLP solver!)
gfcn_nlp_mx.init()
gfcn_nlp = gfcn_nlp_mx.expand([V_sx])
g_sx = gfcn_nlp.outputSX()
  
# Expand objective function as an SX graph (move to NLP solver!)
ffcn_nlp_mx.init()
ffcn_nlp = ffcn_nlp_mx.expand([V_sx])
f_sx = ffcn_nlp.outputSX()
  
# Lagrange multipliers
lam = symbolic("lambda",g_sx.size1())

# Objective function scaling
sigma = symbolic("sigma")

# Lagrangian function (move to NLP solver!)
lfcn = SXFunction([V_sx,lam,sigma], [sigma*f_sx + inner_prod(lam,g_sx)])
  
# Hessian of the Lagrangian
HL = lfcn.hessian()

# Lagrange multipliers
lam_mx = MX("lambda",g_sx.size1())

# Objective function scaling
sigma_mx = MX("sigma")

# Lagrangian function (move to NLP solver!)
lfcn_mx = MXFunction([V,lam_mx,sigma_mx], [sigma_mx*f + inner_prod(lam_mx,gg)])
lfcn_mx.init()

# Gradient of the Lagrangian
gl  = lfcn_mx.grad()
glfcn = MXFunction([V,lam_mx,sigma_mx],[trans(gl[0])])
glfcn.init()

# Hessian of the Lagrangian
HL_mx = glfcn.jacobian()
HL_mx.init()

#raise(Exception("oj"))

## ----
## SOLVE THE NLP
## ----
  
# Allocate an NLP solver
solver = IpoptSolver(ffcn_nlp_mx,gfcn_nlp_mx,HL_mx)

# Set options
solver.setOption("tol",1e-10)
#solver.setOption("hessian_approximation","limited-memory")
#solver.setOption("derivative_test","first-order")
solver.setOption("max_iter",1000)

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,NLP_X_INIT)

# Bounds on x
solver.setInput(vars_lb,NLP_LBX)
solver.setInput(vars_ub,NLP_UBX)

# Bounds on g
solver.setInput(NP.concatenate(lbg),NLP_LBG)
solver.setInput(NP.concatenate(ubg),NLP_UBG)

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", float(solver.output(NLP_COST))

# Retrieve the solution
v_opt = NP.array(solver.output(NLP_X_OPT))

# Get values at the beginning of each finite element
xx_opt = v_opt[0::(K+1)*nx+nu]
yy_opt = v_opt[1::(K+1)*nx+nu]
ll_opt = v_opt[2::(K+1)*nx+nu]
uu_opt = v_opt[(K+1)*nx::(K+1)*nx+nu]

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(uu_opt)
plt.title("u_opt")

plt.figure(2)
plt.clf()
plt.plot(xx_opt)
plt.title("x_opt")

plt.figure(3)
plt.clf()
plt.plot(yy_opt)
plt.title("y_opt")

plt.figure(4)
plt.clf()
plt.plot(ll_opt)
plt.title("l_opt")

plt.show()

