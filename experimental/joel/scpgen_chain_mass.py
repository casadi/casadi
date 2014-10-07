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
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

""" 
Optimal control of a chain mass model
L. Wirsching, H.G. Bock, M. Diehl, Fast NMPC of a chain of masses connected by springs
Computer Aided Control System Design, 2006 IEEE International Conference on Control Applications, 2006

@author Joel Andersson, Milan Vukov, K.U. Leuven 2013
"""

# Start of the chain
pStart = [0., 0., 0.]

# End of the chain in steady state
pEnd = [1., 0., 0.]

# Model parameters
g = 9.81           # Gravitational acceleration [N/kg]
L = 0.033          # Rest length of the springs [m]
D = 1.0            # Spring constant [N/m] (0.1 in Wirsching2006)
m = 0.03           # Ball mass [kg]

# Number of balls
N = 9

# Position and velocities of the masses
p = ssym("p", 3,N+1)
pdot = ssym("pdot", 3,N)  

# Initial guess (uniformly distributed and in rest)
pInit = DMatrix.zeros(3,N+1)
pInit[2,:] = linspace(0,1.0,N+1).T
pdotInit = DMatrix.zeros(3,N)

# Control vector (velocity of the last ball)
u = ssym("u",3)

# Forces between all the balls
F = SX.zeros(3,N+1)
for i in range(N+1):
  # Spring vector
  if i==0:
    d = p[:,0] - pStart
  else:
    d = p[:,i] - p[:,i-1]

  # Length of the spring (add small constant to avoid square root of zero)
  L_d = sqrt(inner_prod(d,d) + 1e-10)

  # Force acting on ball i+1 
  F[:,i] = D*(1.0 - L/L_d) * d

# Acceleration of the masses
a = (F[:,1:]-F[:,0:-1])/m

# Add gravity
a[2,:] -= g

# State vector with initial guess as well as differential equations
x = SX.zeros(0,1)
xInit = DMatrix.zeros(0,1)
f = SX.zeros(0,1)
for i in range(N):
  # Position
  f.append(pdot[:,i])
  x.append(p[:,i])
  xInit.append(pInit[:,i])

  # Velocity
  f.append(a[:,i])
  x.append(pdot[:,i])
  xInit.append(pdotInit[:,i])

# Last ball is controlled
f.append(u)
x.append(p[:,N])
xInit.append(pInit[:,N])

# Define the optimal control problem
nx = x.size()
nu = u.size()

# Weighting factors in the objective
alpha = 25 # 1/s2
beta = 1
gamma = 0.01

# Deviation from the end point
dpEnd = p[:,-1]-pEnd

# Cost function
L = alpha * inner_prod(dpEnd,dpEnd) + beta * inner_prod(vec(p),vec(p)) + gamma * inner_prod(u,u)

# Number of shooting intervals
nk = 20

# Time horizon
T = 8.0

# Control bounds
u_max = 1e-3 # 1/s

# ODE function
ffcn = SXFunction(daeIn(x=x,p=u),daeOut(ode=f,quad=L))

# Create an integrator
Ffcn = CVodesIntegrator(ffcn)
Ffcn.setOption("abstol",1e-10) # tolerance
Ffcn.setOption("reltol",1e-10) # tolerance
Ffcn.setOption("tf",T/nk) # final time
Ffcn.init()

# All controls
U = msym("U",3*nk)

# The initial state
X  = xInit

# Const function
J = 0

# Build a graph of integrator calls
for k in range(nk):
  X,Q,_,_ = Ffcn.call(integratorIn(x0=X,p=U[3*k:3*(k+1)]))
  J += Q
  
# Objective function
F = MXFunction([U],[J])

# Allocate an NLP solver
solver = IpoptSolver(F)
solver.init()

# Set bounds and initial guess
solver.setInput(-u_max, "lbx")
solver.setInput( u_max, "ubx")
solver.setInput( 0.,    "x0")

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = NP.array(solver.getOutput("x"))

