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
from numpy import *
import matplotlib.pyplot as plt

# Time 
t = SX("t")

# Differential states
s = SX("s"); v = SX("v"); m = SX("m")
y = [s,v,m]

# Control
u = SX("u")

alpha = 0.05 # friction
beta = 0.1 # fuel consumption rate

# Differential equation
sdot = v
vdot = (u-alpha*v*v)/m
mdot = -beta*u*u
rhs = [sdot,vdot,mdot]

# ODE right hand side
ffcn = SXFunction([[t],y,[u]],[rhs])
ffcn.setOption("name","ODE right hand side")

# Explicit integrator (CVODES)
integrator = CVodesIntegrator(ffcn)

# Set options
integrator.setOption("fsens_err_con",True)
integrator.setOption("quad_err_con",True)
integrator.setOption("abstol",1e-6)
integrator.setOption("reltol",1e-6)

# Initialize the integrator
integrator.init()

# Time horizon
T = 10.0

# Shooting length
nu = 100 # Number of control segments
DT = T/nu

# Initial position, speed and mass
s0 = 0 # initial position
v0 = 0 # initial speed
m0 = 1 # initial mass

# control for all segments
U = MX("U",nu)

# Dummy input corresponding to the state derivative
xdot = MX([0,0,0]) 

# Integrate over all intervals
X=MX([s0,v0,m0])
T0 = MX(0) # Beginning of time interval (changed from k*DT due to probable Sundials bug)
TF = MX(DT) # End of time interval (changed from (k+1)*DT due to probable Sundials bug)
for k in range(nu):
  # build up a graph with function calls
  X = integrator([T0,TF,X,U[k],xdot])  

# Objective function
F = inner_prod(U,U)

# Terminal constraints
G = vertcat((X[0],X[1]))

# Create the NLP
ffcn = MXFunction([U],[F]) # objective function
gfcn = MXFunction([U],[G]) # constraint function

# Allocate an NLP solver
solver = IpoptSolver(ffcn,gfcn)

# Set options
solver.setOption("tol",1e-10)
solver.setOption("hessian_approximation","limited-memory");

# initialize the solver
solver.init()

# Bounds on u and initial condition
Umin = nu * [-10] # lower bound
solver.setInput(Umin,"lbx")

Umax = nu * [10]  # upper bound
solver.setInput(Umax,"ubx")

Usol = nu * [0.4] # initial guess
solver.setInput(Usol,"x0")

# Bounds on g
Gmin = Gmax = [10, 0]
solver.setInput(Gmin,"lbg")
solver.setInput(Gmax,"ubg")

# Solve the problem
solver.solve()

# Get the solution
uopt = solver.getOutput("x")

# Plot the optimal trajectory
tgrid = linspace(0,T,nu+1)
tgrid_u = linspace(0,T,nu)
plt.figure(1)
plt.clf()
plt.ylabel('Optimal control')
plt.xlabel('time')
plt.plot(tgrid_u,uopt) 

x = [0, 0, 1]
sopt = [x[0]]
vopt = [x[1]]
mopt = [x[2]]
for k in range(nu):
  integrator.setInput(k*DT,"t0")
  integrator.setInput((k+1)*DT,"tf")
  integrator.setInput(uopt[k],"p")
  integrator.setInput(x,"x0")
  integrator.evaluate()
  x = integrator.getOutput()
  sopt.append(x[0])
  vopt.append(x[1])
  mopt.append(x[2])

plt.figure(2)
plt.clf()
plt.subplot(3,1,1)
plt.ylabel('Distance')
plt.xlabel('time')
plt.plot(tgrid,sopt) 

plt.subplot(3,1,2)
plt.ylabel('Velocity')
plt.xlabel('time')
plt.plot(tgrid,vopt) 

plt.subplot(3,1,3)
plt.ylabel('Mass')
plt.xlabel('time')
plt.plot(tgrid,mopt) 

plt.show()
  
