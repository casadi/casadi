#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
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
from casadi.tools import *

import numpy as NP
import matplotlib.pyplot as plt

nk = 20    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
t  = SX.sym("t")    # time
u  = SX.sym("u")    # control

states = struct_symSX([
            entry('x',shape=2),    #  vdp oscillator states
            entry('L')             #  helper state: Langrange integrand
         ])

# Create a structure for the right hand side
rhs = struct_SX(states)
x = states['x']
rhs["x"] = vertcat((1 - x[1]*x[1])*x[0] - x[1] + u, x[0])
rhs["L"] = x[0]*x[0] + x[1]*x[1] + u*u

# ODE right hand side function
f = Function('f', [t,states,u],[rhs])

# Objective function (meyer term)
m = Function('m', [t,states,u],[states["L"]])

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

# Dimensions
nx = 3
nu = 1

# Choose collocation points
tau_root = [0] + collocation_points(3,"radau")

# Degree of interpolating polynomial
d = len(tau_root)-1

# Size of the finite elements
h = tf/nk

# Coefficients of the collocation equation
C = NP.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = NP.zeros(d+1)

# Dimensionless time inside one control interval
tau = SX.sym("tau")

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

  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn = Function('lfcn', [tau],[L])
  D[j] = lfcn(1.0)

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  tfcn = Function('tfcn', [tau],[tangent(L,tau)])
  for r in range(d+1):
    C[j,r] = tfcn(tau_root[r])

# Structure holding NLP variables
V = struct_symMX([
      (
       entry("X",repeat=[nk+1,d+1],struct=states),
       entry("U",repeat=[nk],shape=nu)
      )
    ])

vars_lb   = V()
vars_ub   = V()
vars_init = V()

# Set states and its bounds
vars_init["X",:,:] = repeated(repeated(x_init))
vars_lb["X",:,:]   = repeated(repeated(x_min))
vars_ub["X",:,:]   = repeated(repeated(x_max))

# Set controls and its bounds
vars_init["U",:] = repeated(u_init)
vars_lb["U",:]   = repeated(u_min)
vars_ub["U",:]   = repeated(u_max)

# State at initial time
vars_lb["X",0,0] = xi_min
vars_ub["X",0,0] = xi_max

# State at end time
vars_lb["X",-1,0] = xf_min
vars_ub["X",-1,0] = xf_max

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
      xp_jk += C[r,j]*V["X",k,r]

    # Add collocation equations to the NLP
    fk = f(T[k][j], V["X",k,j], V["U",k])
    g.append(h*fk - xp_jk)
    lbg.append(NP.zeros(nx)) # equality constraints
    ubg.append(NP.zeros(nx)) # equality constraints

  # Get an expression for the state at the end of the finite element
  xf_k = 0
  for r in range(d+1):
    xf_k += D[r]*V["X",k,r]

  # Add continuity equation to NLP
  g.append(V["X",k+1,0] - xf_k)
  lbg.append(NP.zeros(nx))
  ubg.append(NP.zeros(nx))

# Concatenate constraints
g = vertcat(*g)

# Objective function
f = m(T[nk-1][d],V["X",nk,0],V["U",nk-1])

# NLP
nlp = {'x':V, 'f':f, 'g':g}

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
opt = V(res["x"])

# Get values at the beginning of each finite element
x0_opt = opt["X",:,0,"x",0]
x1_opt = opt["X",:,0,"x",1]
x2_opt = opt["X",:,0,"L"]

u_opt = opt["U",:,0]

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
