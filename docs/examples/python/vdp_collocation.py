#
#     MIT No Attribution
#
#     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# -*- coding: utf-8 -*-
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

# Degree of interpolating polynomial
d = 3

# Choose collocation points
tau_root = [0] + collocation_points(d, "radau")

# Coefficients of the collocation equation
C = NP.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = NP.zeros(d+1)

# Coefficients of the quadrature function
F = NP.zeros(d+1)

# Construct polynomial basis
for j in range(d+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  p = NP.poly1d([1])
  for r in range(d+1):
    if r != j:
      p *= NP.poly1d([1, -tau_root[r]]) / (tau_root[j]-tau_root[r])
  
  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D[j] = p(1.0)

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = NP.polyder(p)
  for r in range(d+1):
    C[j,r] = pder(tau_root[r])

  # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = NP.polyint(p)
  F[j] = pint(1.0)

# Control discretization
nk = 20

# End time
tf = 10.0  

# Size of the finite elements
h = tf/nk

# All collocation time points
T = NP.zeros((nk,d+1))
for k in range(nk):
  for j in range(d+1):
    T[k,j] = h*(k + tau_root[j])

# Declare variables (use scalar graph)
t  = SX.sym("t")    # time
u  = SX.sym("u")    # control
x  = SX.sym("x",2)  # state

# ODE rhs function and quadratures
xdot = vertcat((1 - x[1]*x[1])*x[0] - x[1] + u, \
               x[0])
qdot = x[0]*x[0] + x[1]*x[1] + u*u
f = Function('f', [t,x,u],[xdot, qdot])

# Control bounds
u_min = -0.75
u_max = 1.0
u_init = 0.0

u_lb = NP.array([u_min])
u_ub = NP.array([u_max])
u_init = NP.array([u_init])

# State bounds and initial guess
x_min =  [-inf, -inf]
x_max =  [ inf,  inf]
xi_min = [ 0.0,  1.0]
xi_max = [ 0.0,  1.0]
xf_min = [ 0.0,  0.0]
xf_max = [ 0.0,  0.0]
x_init = [ 0.0,  0.0]

# Dimensions
nx = 2
nu = 1

# Total number of variables
NX = nk*(d+1)*nx      # Collocated states
NU = nk*nu              # Parametrized controls
NXF = nx                # Final state
NV = NX+NU+NXF

# NLP variable vector
V = MX.sym("V",NV)
  
# All variables with bounds and initial guess
vars_lb = NP.zeros(NV)
vars_ub = NP.zeros(NV)
vars_init = NP.zeros(NV)
offset = 0

# Get collocated states and parametrized control
X = NP.resize(NP.array([],dtype=MX),(nk+1,d+1))
U = NP.resize(NP.array([],dtype=MX),nk)
for k in range(nk):  
  # Collocated states
  for j in range(d+1):
    # Get the expression for the state vector
    X[k,j] = V[offset:offset+nx]
    
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
X[nk,0] = V[offset:offset+nx]
vars_lb[offset:offset+nx] = xf_min
vars_ub[offset:offset+nx] = xf_max
vars_init[offset:offset+nx] = x_init
offset += nx
  
# Constraint function for the NLP
g = []
lbg = []
ubg = []

# Objective function
J = 0

# For all finite elements
for k in range(nk):
  
  # For all collocation points
  for j in range(1,d+1):
        
    # Get an expression for the state derivative at the collocation point
    xp_jk = 0
    for r in range (d+1):
      xp_jk += C[r,j]*X[k,r]
      
    # Add collocation equations to the NLP
    fk,qk = f(T[k,j], X[k,j], U[k])
    g.append(h*fk - xp_jk)
    lbg.append(NP.zeros(nx)) # equality constraints
    ubg.append(NP.zeros(nx)) # equality constraints

    # Add contribution to objective
    J += F[j]*qk*h

  # Get an expression for the state at the end of the finite element
  xf_k = 0
  for r in range(d+1):
    xf_k += D[r]*X[k,r]

  # Add continuity equation to NLP
  g.append(X[k+1,0] - xf_k)
  lbg.append(NP.zeros(nx))
  ubg.append(NP.zeros(nx))
  
# Concatenate constraints
g = vertcat(*g)
  
# NLP
nlp = {'x':V, 'f':J, 'g':g}

## ----
## SOLVE THE NLP
## ----

# Set options
opts = {}
opts["expand"] = True
#opts["ipopt.max_iter"] = 4
opts["ipopt.linear_solver"] = 'ma27'

# Allocate an NLP solver
solver = nlpsol("solver", "ipopt", nlp, opts)
arg = {}
  
# Initial condition
arg["x0"] = vars_init

# Bounds on x
arg["lbx"] = vars_lb
arg["ubx"] = vars_ub

# Bounds on g
arg["lbg"] = NP.concatenate(lbg)
arg["ubg"] = NP.concatenate(ubg)

# Solve the problem
res = solver(**arg)

# Print the optimal cost
print("optimal cost: ", float(res["f"]))

# Retrieve the solution
v_opt = NP.array(res["x"])

# Get values at the beginning of each finite element
x0_opt = v_opt[0::(d+1)*nx+nu]
x1_opt = v_opt[1::(d+1)*nx+nu]
u_opt = v_opt[(d+1)*nx::(d+1)*nx+nu]
tgrid = NP.linspace(0,tf,nk+1)
tgrid_u = NP.linspace(0,tf,nk)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid,x0_opt,'--')
plt.plot(tgrid,x1_opt,'-.')
plt.step(tgrid_u,u_opt,'-')
plt.title("Van der Pol optimization")
plt.xlabel('time')
plt.legend(['x[0] trajectory','x[1] trajectory','u trajectory'])
plt.grid()
plt.show()

