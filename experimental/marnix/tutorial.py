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
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt

# CasADi
from casadi import *

# Variables in the state equation
x = ssym("x",2)
u = ssym("u")
w = ssym("w")

# Dimensions
nx = 2
nu = 1
nw = 1
nz = 1

# Regularization
x_reg = 1.0
u_reg = 0.01

# State equation
f = [ 1.0*x[0] + 0.1*x[1],\
     -0.5*x[0] + 0.9*x[1] + 0.1*u]
f = vertcat(f)

ffcn = SXFunction([x,u,w],[f])
ffcn.init()

# Output equation
hz = x[0]**2 + x[1]**2

hzfcn = SXFunction([x,u,w],[hz])
hzfcn.init()

# Time discretization
N = 1000

# State
X = ssym("X",2,N)

# Control
U = ssym("U",1,N)

# Outputs to be constrained
Z = ssym("Z",1,N)

# Assemble the NLP variable vector and bounds
v = []
v_min = []
v_max = []
v_init = []
for k in range(N):
  v += [X[:,k]]
  v += [U[:,k]]
  v += [Z[:,k]]

  v_min += [-inf*ones(2)] # X
  v_min += [-inf*ones(1)] # U
  v_min += [-inf*ones(1)] # Z

  v_max += [ inf*ones(2)]
  v_max += [ inf*ones(1)]
  v_max += [ 1.0*ones(1)]

  v_init += [ 0.0*ones(2)]
  v_init += [ 0.0*ones(1)]
  v_init += [ 0.0*ones(1)]

v = vertcat(v)
v_min = NP.concatenate(v_min)
v_max = NP.concatenate(v_max)
v_init = NP.concatenate(v_init)

# Objective function weighting
Q = concatenate((x_reg*ones(N),u_reg*ones(N)))
Q = diag(Q)

# Outputs to be minimized
W = ones(N)

# Residual
y_bar = W-trans(X[1,:])
y_bar = vertcat((y_bar,trans(U)))

# Objective
J = dot(trans(y_bar),dot(Q,y_bar))

# Objective function
Jfcn = SXFunction([v],[J])

# Initial state
x_init = NP.array([0,0])

# Constraint function
g = []
g += [X[:,0]-x_init]
for k in range(0,N-1):
  # Evaluate the state equation symbolically
  [x_next] = ffcn.eval([X[:,k],U[:,k],W[k]])

  # Append to the constraint expression
  g += [X[:,k+1]-x_next]

  # Evaluate the constraint function symbolically
  [z_k] = hzfcn.eval([X[:,k],U[:,k],W[k]])
  
  # Append to the constraint expression
  g += [Z[:,k]-z_k]

# Add the final point
# Evaluate the constraint function symbolically
[z_final] = hzfcn.eval([X[:,N-1],U[:,N-1],W[N-1]])

# Append to the constraint expression
g += [Z[:,N-1]-z_final]

# Concatenate
g = vertcat(g)
  
# Constraint function
gfcn = SXFunction([v],[g])
  
# Allocate an IPOPT instance
solver = IpoptSolver(Jfcn,gfcn)

# Set IPOPT options
solver.setOption("tol",1e-10)
solver.setOption("generate_hessian",True)
#solver.setOption("hessian_approximation","limited-memory")

# Initialize IPOPT
solver.init()

# Pass bounds
solver.setInput(v_min,"lbx") # lower variable bounds
solver.setInput(v_max,"ubx") # upper variable bounds
solver.setInput(v_init,"x0") # variable initial guess
solver.setInput(NP.zeros(g.size()),"lbg") # equality constraints
solver.setInput(NP.zeros(g.size()),"ubg") # equality constraints

# Solve the NLP
solver.solve()

# Get the solution
v_opt = NP.array(solver.getOutput("x"))
x_opt = zeros(X.shape)
u_opt = zeros(U.shape)
z_opt = zeros(Z.shape)
vk = 0
for k in range(N):
  # Get the state
  for i in range(nx):
    x_opt[i,k] = v_opt[vk+i]
  vk += nx
  
  # Get the inputs
  for i in range(nu):
    u_opt[i,k] = v_opt[vk+i]
  vk += nu
  
  # Skip the constrained output
  for i in range(nz):
    z_opt[i,k] = v_opt[vk+i]
  vk += nz

# Time grid
t_grid = NP.linspace(0,1,N)

# Plot the solution
plt.figure(1)
plt.clf()
plt.plot(x_opt[0,:],x_opt[1,:],'*')
plt.xlabel('x_0')
plt.ylabel('x_1')

# also plot a circle
phi = NP.linspace(0,2*pi,100)
plt.plot(sin(phi),cos(phi),'-')
plt.axis('equal')

plt.figure(2)
plt.clf()
plt.plot(t_grid,u_opt.T)
plt.xlabel('time')
plt.ylabel('u')

plt.figure(3)
plt.clf()
plt.plot(t_grid,x_opt[0,:].T)
plt.xlabel('time')
plt.ylabel('x0')

plt.figure(4)
plt.clf()
plt.plot(t_grid,x_opt[1,:].T)
plt.xlabel('time')
plt.ylabel('x1')

plt.show()
