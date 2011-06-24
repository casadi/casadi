# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt

# Excercise 2, chapter 10 from Larry Biegler's book 
# Joel Andersson, K.U. Leuven 2010

print "program started"

# Time 
t = SX("t")
  
# Differential states
Ls = SX("Ls")    # mean crystal size
Nc = SX("Nc")    # number of nuclei per liter of solvent
L = SX("L")      # total length of crystals per liter of solvent
Ac = SX("Ac")    # total surface area of the crystals per liter of solvent
Vc = SX("Vc")    # total volume of the crysals per liter of solvent
Mc = SX("Mc")    # total mass of the crystals
Cc = SX("Cc")    # solute concentration
Tc = SX("Tc")    # cystillizer temperature

# State vector
x = [Ls, Nc, L, Ac, Vc, Mc, Cc, Tc]
  
# Bounds on the states
x_lb = NP.repeat(-inf,len(x))
x_ub = NP.repeat( inf,len(x))

# Initial values
Ls_init = 0.0005
Nc_init = 0
L_init = 0
Ac_init = 0
Vc_init = 0
Mc_init = 2.0
Cc_init = 5.4
Tc_init = 75
x_init = NP.array([Ls_init,Nc_init,L_init,Ac_init,Vc_init,Mc_init,Cc_init,Tc_init])

# Controls
Tj = SX("Tj") # jacket temperature
Tj_lb = 10;  Tj_ub = 60;  Tj_init = 50

# Constants
Vs = 300. # volume of the solvent
W = 2025. # the total mass in the crysallizer
a = [-66.4309, 2.8604, -0.022579, 6.7117e-5]
b = [16.08852, -2.708263, 0.0670694, -3.5685e-4]
Kg = 0.00418
Bn = 385.
Cp = 0.4
Kc = 35.
Ke = 377.
eta1 = 1.1
eta2 = 5.72
Ls0 = 5e-4 # initial crystal size
L0 = 5e-5 # nucleate crystal size
Ws0 = 2. # weight of the seed crystals
rho = 1.58 # specific gravity of crystals
alpha = 0.2 # shape factor for area of crystals
beta = 1.2 # shape factor for volume of crystals

# Time horizon
tf = 25.

# Dependent variables
C_bar = 100.*Cc/(1.35+Cc)
Tequ = a[0] + a[1]*C_bar + a[2]*C_bar*C_bar + a[3]*C_bar*C_bar*C_bar # equilibrium temperature
Ta = b[0] + b[1]*C_bar + b[2]*C_bar*C_bar + b[3]*C_bar*C_bar*C_bar # lower bound of Tj

# degree of supercooling:
# DeltaT = fmax(0,Tequ-Tc)   # Original formulation
# DeltaT = fmax(1e-8,Tequ-Tc) # "epsilon" to avoid divide by zero
# DeltaT = log(exp(1e-8)+exp(Tequ-Tc))   # Log sum exp
epsilon = 1e-3
DD = Tequ-Tc
DeltaT = (sqrt( DD*DD + epsilon*epsilon )  + DD ) / 2

# Differential equations
Ls_dot = Kg * sqrt(Ls) * DeltaT**eta1
Nc_dot = Bn * DeltaT**eta2
L_dot = Nc * Ls_dot + L0 * Nc_dot
Ac_dot = 2 * alpha * Nc * Ls_dot + L0**2 * Nc_dot
Vc_dot = 3 * beta * Ac * Ls_dot + L0**3 * Nc_dot
Mc_dot = 3 * Ws0/Ls0**3 * Ls**2 * Ls_dot + rho*Vs*Vc_dot
Cc_dot = -1 / Vs * Mc_dot
Tc_dot = (Kc*Mc_dot - Ke*(Tc-Tj))/(W*Cp)

# Control vector
u = [Tj]
u_init = NP.array([Tj_init])
u_lb = NP.array([Tj_init])
u_ub = NP.array([Tj_init])
#u_lb = [Tj_lb]
#u_ub = [Tj_ub]

# ODE
xdot = [Ls_dot, Nc_dot, L_dot, Ac_dot, Vc_dot, Mc_dot, Cc_dot, Tc_dot]

# Right hand side of the ODE
ffcn = SXFunction([[t],x,u],[xdot])
ffcn.init()
  
# Objective function (meyer term)
mfcn = SXFunction([[t],x,u],[[-Ls]])
mfcn.init()

# Nonlinear constraint function
cfcn = SXFunction([[t],x,u],[[Tj-Ta]])
cfcn.init()

# Dimensions
nx = len(x)
nu = len(u)


print "modelling done"


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
N = 30

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
V = symbolic("V",NV)
V2 = MX("V",NV)

# Collocated states
X = []
X2 = []

# Collocated control (piecewice constant)
U = []
U2 = []
  
# All variables with bounds and initial guess
vars_lb = NP.zeros(NV)
vars_ub = NP.zeros(NV)
vars_init = NP.zeros(NV)
vi = 0

# Loop over the finite elements
for i in range(N):
  # Parametrized controls
  U.append(V[vi:vi+nu])
  U2.append(V2[vi:vi+nu])
  vars_lb[vi:vi+nu] = u_lb
  vars_ub[vi:vi+nu] = u_ub
  vars_init[vi:vi+nu] = u_init
  vi += nu
  
  # Collocated states
  Xi = []
  Xi2 = []
  for j in range(K+1):
    Xi.append(V[vi:vi+nx])
    Xi2.append(V2[vi:vi+nx])
    vars_init[vi:vi+nx] = x_init
    if i==0 and j==0:
      vars_lb[vi:vi+nx] = x_init
      vars_ub[vi:vi+nx] = x_init
    else:
      vars_lb[vi:vi+nx] = x_lb
      vars_ub[vi:vi+nx] = x_ub
    vi += nx
        
  X.append(Xi)
  X2.append(Xi2)
  
# State at end time
XF = V[vi:vi+nx]
XF2 = V2[vi:vi+nx]
vars_lb[vi:vi+nx] = x_lb
vars_ub[vi:vi+nx] = x_ub
vars_init[vi:vi+nx] = x_init
vi += nx
  
# Constraint function for the NLP
g = []
g2 = []
lbg = []
ubg = []

for i in range(N):
  for k in range(1,K+1):
    # augmented state vector
    y_ik = [[SX(T[i][k])], X[i][k],  U[i]]
    y_ik2 = [MX(T[i][k]), X2[i][k], U2[i]]
    
    # Add collocation equations to NLP
    [temp] = ffcn.eval(y_ik)
    [temp2] = ffcn.call(y_ik2)
    for j in range(temp.size()):
      temp[j] *= h
      temp2[j] *= h
    for j in range (K+1):
      for l in range(temp.size()):
        temp[l] -= X[i][j][l]*C[j][k]
        temp2[l] -= X2[i][j][l]*C[j][k]
      
    g += [temp]
    g2 += [temp2]
    lbg += nx * [0.] # equality constraints
    ubg += nx * [0.] # equality constraints
      
    # Add nonlinear constraints
    [temp] = cfcn.eval(y_ik)
    [temp2] = cfcn.call(y_ik2)
    g += [temp]
    g2 += [temp2]
    lbg += [0.] # correct!
    ubg += [inf]

  # Add continuity equation to NLP
  if i<N-1:
    temp = SXMatrix(X[i+1][0])
    temp2 = MX(X2[i+1][0])
  else:
    temp = SXMatrix(XF)
    temp2 = MX(XF2)
    
  for j in range(K+1):
    for l in range(nx):
      temp[l] -= D[j]*X[i][j][l]
      temp2[l] -= D[j]*X2[i][j][l]

  g += [temp]
  g2 += [temp2]
  lbg += nx*[0.]
  ubg += nx*[0.]
  
# Nonlinear constraint function
g_all = vertcat(g)
g_all2 = vertcat(g2)

gfcn_nlp = SXFunction([V],[g_all])
gfcn_nlp2 = MXFunction([V2],[g_all2])
gfcn_nlp2.init()
gfcn_nlp3 = SXFunction(gfcn_nlp2)
  
# Objective function of the NLP
y_f = [[SX(T[N-1][K])],XF,U[N-1]]
f = mfcn.eval(y_f)[0][0]
ffcn_nlp = SXFunction([V], [[f]])
  
# Hessian of the Lagrangian:
# Lagrange multipliers
lam = create_symbolic("lambda",g_all.size1())

# Objective function scaling
sigma = SX("sigma")

# Lagrangian function
lfcn_input = [V,lam,[sigma]]
lfcn_output = sigma*f
for j in range(len(g)):
  lfcn_output += lam[j]*g_all.at(j)
lfcn = SXFunction(lfcn_input, [[lfcn_output]])
  
# Hessian of the Lagrangian
HL = lfcn.hessian()

## ----
## SOLVE THE NLP
## ----
  
# Allocate an NLP solver
solver = IpoptSolver(ffcn_nlp,gfcn_nlp3,HL)

# Set options
solver.setOption("tol",1e-6)
#solver.setOption("hessian_approximation","limited-memory")
#solver.setOption("derivative_test","first-order")
solver.setOption("max_iter",50)

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,NLP_X_INIT)

# Bounds on x
solver.setInput(vars_lb,NLP_LBX)
solver.setInput(vars_ub,NLP_UBX)

# Bounds on g
solver.setInput(lbg,NLP_LBG)
solver.setInput(ubg,NLP_UBG)

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", solver.output(NLP_COST)[0]

# Get the solution
vars_sol = solver.output(NLP_X_OPT)

## ----
## SAVE SOLUTION TO DISK
## ----

# Get the optimal solution
Tj_opt = N * [0.]
Ls_opt = (N*(K+1)) * [0.]
Nc_opt = (N*(K+1)) * [0.]
L_opt = (N*(K+1)) * [0.]
Ac_opt = (N*(K+1)) * [0.]
Vc_opt = (N*(K+1)) * [0.]
Mc_opt = (N*(K+1)) * [0.]
Cc_opt = (N*(K+1)) * [0.]
Tc_opt = (N*(K+1)) * [0.]
t_opt = (N*(K+1)) * [0.]
ind = 0 # index of var_sol
for i in range(N):
  Tj_opt[i] = vars_sol[ind]; ind += 1
  for j in range(K+1):
    ij = (K+1)*i+j
    Ls_opt[ij] = vars_sol[ind]; ind += 1
    Nc_opt[ij] = vars_sol[ind]; ind += 1
    L_opt[ij] = vars_sol[ind]; ind += 1
    Ac_opt[ij] = vars_sol[ind]; ind += 1
    Vc_opt[ij] = vars_sol[ind]; ind += 1
    Mc_opt[ij] = vars_sol[ind]; ind += 1
    Cc_opt[ij] = vars_sol[ind]; ind += 1
    Tc_opt[ij] = vars_sol[ind]; ind += 1
    t_opt[ij] = T[i][j]
  
# plot to screen
plt.figure(1)
plt.clf()
plt.plot(t_opt,Ls_opt)
plt.xlabel("Time (h)")
plt.ylabel("Mean Chrystal Size (m)")

plt.figure(2)
plt.clf()
plt.plot(linspace(0,tf,N),Tj_opt)
plt.xlabel("Time (h)")
plt.ylabel("Jacket Temperature (C)")

#plt.figure(3)
#plt.clf()
#plt.plot(t_opt,Ls_opt)


# show the plots
plt.show()
