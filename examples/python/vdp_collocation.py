#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
# 
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
# 
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# 
# 
# -*- coding: utf-8 -*-
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

nk = 20    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
t  = ssym("t")    # time
u  = ssym("u")    # control
x  = ssym("x",3)  # state

# ODE right hand side function
rhs = vertcat([(1 - x[1]*x[1])*x[0] - x[1] + u, \
               x[0], \
               x[0]*x[0] + x[1]*x[1] + u*u])
f = SXFunction([t,x,u],[rhs])

# Objective function (meyer term)
m = SXFunction([t,x,u],[x[2]])

# Control bounds
u_min = -0.75
u_max = 1.0
u_init = 0.0

u_lb = NP.array([u_min])
u_ub = NP.array([u_max])
u_init = NP.array([u_init])

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

# Choose collocation points
tau_root = collocationPoints(3,"radau")

# Degree of interpolating polynomial
d = len(tau_root)-1

# Size of the finite elements
h = tf/nk

# Coefficients of the collocation equation
C = NP.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = NP.zeros(d+1)

# Dimensionless time inside one control interval
tau = ssym("tau")
  
# All collocation time points
T = NP.zeros((nk,d+1))
for k in range(nk):
  for j in range(d+1):
    T[k,j] = h*(k + tau_root[j])

# For all collocation points
for j in range(d+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  L = 1
  for r in range(d+1):
    if r != j:
      L *= (tau-tau_root[r])/(tau_root[j]-tau_root[r])
  lfcn = SXFunction([tau],[L])
  lfcn.init()
  
  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn.setInput(1.0)
  lfcn.evaluate()
  D[j] = lfcn.output()

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  for r in range(d+1):
    lfcn.setInput(tau_root[r])
    lfcn.setFwdSeed(1.0)
    lfcn.evaluate(1,0)
    C[j,r] = lfcn.fwdSens()

# Total number of variables
NX = nk*(d+1)*nx      # Collocated states
NU = nk*nu              # Parametrized controls
NXF = nx                # Final state
NV = NX+NU+NXF

# NLP variable vector
V = MX("V",NV)
  
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

# For all finite elements
for k in range(nk):
  
  # For all collocation points
  for j in range(1,d+1):
        
    # Get an expression for the state derivative at the collocation point
    xp_jk = 0
    for r in range (d+1):
      xp_jk += C[r,j]*X[k,r]
      
    # Add collocation equations to the NLP
    [fk] = f.call([T[k,j], X[k,j], U[k]])
    g.append(h*fk - xp_jk)
    lbg.append(NP.zeros(nx)) # equality constraints
    ubg.append(NP.zeros(nx)) # equality constraints

  # Get an expression for the state at the end of the finite element
  xf_k = 0
  for r in range(d+1):
    xf_k += D[r]*X[k,r]

  # Add continuity equation to NLP
  g.append(X[k+1,0] - xf_k)
  lbg.append(NP.zeros(nx))
  ubg.append(NP.zeros(nx))
  
# Concatenate constraints
g = vertcat(g)

# Objective function
[f] = m.call([T[nk-1,d],X[nk,0],U[nk-1]])
  
# NLP
nlp = MXFunction(nlpIn(x=V),nlpOut(f=f,g=g))
  
## ----
## SOLVE THE NLP
## ----
  
# Allocate an NLP solver
solver = IpoptSolver(nlp)

# Set options
solver.setOption("expand",True)
#solver.setOption("max_iter",4)

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,"x0")

# Bounds on x
solver.setInput(vars_lb,"lbx")
solver.setInput(vars_ub,"ubx")

# Bounds on g
solver.setInput(NP.concatenate(lbg),"lbg")
solver.setInput(NP.concatenate(ubg),"ubg")

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", float(solver.output("f"))

# Retrieve the solution
v_opt = NP.array(solver.output("x"))

# Get values at the beginning of each finite element
x0_opt = v_opt[0::(d+1)*nx+nu]
x1_opt = v_opt[1::(d+1)*nx+nu]
x2_opt = v_opt[2::(d+1)*nx+nu]
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

