# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import matplotlib.pyplot as plt

nk = 20    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
t  = ssym("t")    # time
u  = ssym("u")    # control
x  = ssym("x",3)  # state

# ODE right hand side function
rhs = [(1 - x[1]*x[1])*x[0] - x[1] + u, \
       x[0], \
       x[0]*x[0] + x[1]*x[1] + u*u]
f = SXFunction([t,x,u],[rhs])

# Objective function (meyer term)
m = SXFunction([t,x,u],[x[2]])

# Control bounds
u_min = -0.75
u_max = 1.0
u_init = 0.0

u_lb = array([u_min])
u_ub = array([u_max])
u_init = array([u_init])

# State bounds and initial guess
x_min =  [-inf, -inf, -inf]
x_max =  [ inf,  inf,  inf]
xi_min = [ 0.0,  1.0,  0.0]
xi_max = [ 0.0,  1.0,  0.0]
xf_min = [ 0.0,  0.0, -inf]
xf_max = [ 0.0,  0.0,  inf]
x_init = [ 0.0,  0.0,  0.0]

# Initialize functions
f.init()  
m.init()

# Dimensions
nx = 3
nu = 1

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
deg = 3

# Radau collocation points
cp = RADAU

# Size of the finite elements
h = tf/nk

# Coefficients of the collocation equation
C = zeros((deg+1,deg+1))

# Coefficients of the continuity equation
D = zeros(deg+1)

# Collocation point
tau = ssym("tau")
  
# All collocation time points
tau_root = collocation_points[cp][deg]
T = zeros((nk,deg+1))
for i in range(nk):
  for j in range(deg+1):
    T[i][j] = h*(i + tau_root[j])

# For all collocation points
for j in range(deg+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  L = 1
  for j2 in range(deg+1):
    if j2 != j:
      L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2])
  lfcn = SXFunction([tau],[L])
  lfcn.init()
  
  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn.setInput(1.0)
  lfcn.evaluate()
  D[j] = lfcn.output()

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  for j2 in range(deg+1):
    lfcn.setInput(tau_root[j2])
    lfcn.setFwdSeed(1.0)
    lfcn.evaluate(1,0)
    C[j][j2] = lfcn.fwdSens()

# Total number of variables
NX = nk*(deg+1)*nx      # Collocated states
NU = nk*nu              # Parametrized controls
NXF = nx                # Final state
NV = NX+NU+NXF

# NLP variable vector
V = MX("V",NV)
  
# All variables with bounds and initial guess
vars_lb = zeros(NV)
vars_ub = zeros(NV)
vars_init = zeros(NV)
offset = 0

# Get collocated states and parametrized control
X = resize(array([],dtype=MX),(nk+1,deg+1))
U = resize(array([],dtype=MX),nk)
for k in range(nk):  
  # Collocated states
  for j in range(deg+1):
    # Get the expression for the state vector
    X[k][j] = V[offset:offset+nx]
    
    # Add the initial condition
    vars_init[offset:offset+nx] = x_init
    
    # Add bounds
    if k==0 and j==0:
      vars_lb[offset:offset+nx] = xi_min
      vars_ub[offset:offset+nx] = xi_max
    else:
      vars_lb[offset:offset+nx] = x_min
      vars_ub[offset:offset+nx] = x_max
    offset += nx
  
  # Parametrized controls
  U[k] = V[offset:offset+nu]
  vars_lb[offset:offset+nu] = u_min
  vars_ub[offset:offset+nu] = u_max
  vars_init[offset:offset+nu] = u_init
  offset += nu
  
# State at end time
X[nk][0] = V[offset:offset+nx]
vars_lb[offset:offset+nx] = xf_min
vars_ub[offset:offset+nx] = xf_max
vars_init[offset:offset+nx] = x_init
offset += nx
  
# Constraint function for the NLP
g = []
lbg = []
ubg = []

# For all finite elements
for k in range(nk):
  
  # For all collocation points
  for j in range(1,deg+1):
        
    # Get an expression for the state derivative at the collocation point
    xp_jk = 0
    for j2 in range (deg+1):
      xp_jk += C[j2][j]*X[k][j2]
      
    # Add collocation equations to the NLP
    [fk] = f.call([T[k][j], X[k][j], U[k]])
    g += [h*fk - xp_jk]
    lbg.append(zeros(nx)) # equality constraints
    ubg.append(zeros(nx)) # equality constraints

  # Get an expression for the state at the end of the finite element
  xf_k = 0
  for j in range(deg+1):
    xf_k += D[j]*X[k][j]

  # Add continuity equation to NLP
  g += [X[k+1][0] - xf_k]
  lbg.append(zeros(nx))
  ubg.append(zeros(nx))
  
# Variable vector (SX)
V_sx = symbolic("V",NV)
  
# Nonlinear constraint function
gfcn_nlp_mx = MXFunction([V],[vertcat(g)])

# Objective function of the NLP
y_f = [MX(T[nk-1][deg]),X[nk][0],U[nk-1]]
[f] = m.call(y_f)
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
lfcn.init()

# Hessian of the Lagrangian
HL = lfcn.hessian()

## ----
## SOLVE THE NLP
## ----
  
# Allocate an NLP solver
solver = IpoptSolver(ffcn_nlp,gfcn_nlp,HL)

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,NLP_X_INIT)

# Bounds on x
solver.setInput(vars_lb,NLP_LBX)
solver.setInput(vars_ub,NLP_UBX)

# Bounds on g
solver.setInput(concatenate(lbg),NLP_LBG)
solver.setInput(concatenate(ubg),NLP_UBG)

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", float(solver.output(NLP_COST))

# Retrieve the solution
v_opt = array(solver.output(NLP_X_OPT))

# Get values at the beginning of each finite element
x0_opt = v_opt[0::(deg+1)*nx+nu]
x1_opt = v_opt[1::(deg+1)*nx+nu]
x2_opt = v_opt[2::(deg+1)*nx+nu]
u_opt = v_opt[(deg+1)*nx::(deg+1)*nx+nu]
tgrid = linspace(0,tf,nk+1)
tgrid_u = linspace(0,tf,nk)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid,x0_opt,'--')
plt.plot(tgrid,x1_opt,'-')
plt.plot(tgrid_u,u_opt,'-.')
plt.title("Van der Pol optimization")
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

