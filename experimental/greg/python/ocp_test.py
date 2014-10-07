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

nk = 20      # Control discretization
tf = 10.0    # End time
coll = True # Use collocation integrator
#coll = False

# Declare variables (use scalar graph)
numStates = 3
numActions = 1

def actionIdx(actionK, timestepK):
  if timestepK > nk - 1:
    raise ValueError("timestepK > nk - 1")
  if timestepK < 0:
    raise ValueError("timestepK < 0")
  if actionK > numActions - 1:
    raise ValueError("actionK > numActions - 1")
  if actionK < 0:
    raise ValueError("actionK < 0")
  return timestepK*numActions + actionK

def stateIdx(stateK, timestepK):
  if timestepK > nk:
    raise ValueError("timestepK > nk")
  if timestepK < 0:
    raise ValueError("timestepK < 0")
  if stateK > numStates - 1:
    raise ValueError("stateK > numStates - 1")
  if stateK < 0:
    raise ValueError("stateK < 0")
  return timestepK*numStates + stateK + numActions*nk

t  = ssym("t")    # time
u  = ssym("u",numActions) # control
x  = ssym("x",numStates)  # state
xp = ssym("xd",numStates) # state derivative

# ODE/DAE residual function
res = [(1 - x[1]*x[1])*x[0] - x[1] + u, \
       x[0], \
       x[0]*x[0] + x[1]*x[1] + u*u] - xp
f = SXFunction([t,x,u,xp],[res])

# Create an integrator
if coll:
  f_d = CollocationIntegrator(f)
  f_d.setOption("number_of_finite_elements",5)
  f_d.setOption("interpolation_order",5)
  f_d.setOption("collocation_scheme","legendre")
  f_d.setOption("implicit_solver",KinsolSolver)
  f_d.setOption("implicit_solver_options",\
    {'linear_solver' : CSparse})
  f_d.setOption("expand_f",True)
else:
  f_d = CVodesIntegrator(f)
  f_d.setOption("abstol",1e-8) # tolerance
  f_d.setOption("reltol",1e-8) # tolerance
  f_d.setOption("steps_per_checkpoint",1000)

f_d.setOption("tf",tf/nk) # final time
f_d.init()

# Total number of variables
nv = numActions*nk + numStates*(nk+1)

# Declare variable vector
V = msym("V", nv)

# Variable bounds initialized to +/- inf
VMIN = -inf*NP.ones(nv)
VMAX = inf*NP.ones(nv)

# Control bounds
for k in range(nk):
  VMIN[actionIdx(0,k)] = -0.75
  VMAX[actionIdx(0,k)] = 1.0

# State bounds
for k in range(nk+1):
  VMIN[stateIdx(0,k)] = -inf
  VMAX[stateIdx(0,k)] = inf

# Initial condition
VMIN[stateIdx(0,0)] = VMAX[stateIdx(0,0)] = 0
VMIN[stateIdx(1,0)] = VMAX[stateIdx(1,0)] = 1
VMIN[stateIdx(2,0)] = VMAX[stateIdx(2,0)] = 0

# Terminal constraint
VMIN[stateIdx(0,nk)] = VMAX[stateIdx(0,nk)] = 0
VMIN[stateIdx(1,nk)] = VMAX[stateIdx(1,nk)] = 0

# Initial solution guess
VINIT = NP.zeros(nv)

# State derivative (only relevant for DAEs)
Xp = msym([0 for k in range(numStates)])

# Constraint function with bounds
g = [];  g_min = []; g_max = []

# Build up a graph of integrator calls
for k in range(nk):
  # Local state/action vectors
  Xk      = vertcat([V[stateIdx( j, k  )] for j in range(numStates)])
  Xk_next = vertcat([V[stateIdx( j, k+1)] for j in range(numStates)])
  Uk      = vertcat([V[actionIdx(j, k  )] for j in range(numActions)])

  # Call the integrator
  [Xk_end,Xp] = f_d.call([Xk,Uk,Xp])
  
  # append continuity constraints
  g.append(Xk_next - Xk_end)
  g_min.append(NP.zeros(Xk.size()))
  g_max.append(NP.zeros(Xk.size()))

# Objective function: L(T)
F = MXFunction([V],[V[stateIdx(numStates-1,nk)]])

# Terminal constraints: 0<=[x(T);y(T)]<=0
G = MXFunction([V],[vertcat(g)])

# Create NLP solver instance
solver = IpoptSolver(F,G)

#solver.setOption("verbose",True)
solver.init()

# Set bounds and initial guess
solver.setInput(VMIN,  "lbx")
solver.setInput(VMAX,  "ubx")
solver.setInput(VINIT, "x0")
solver.setInput(NP.concatenate(g_min),"lbg")
solver.setInput(NP.concatenate(g_max),"ubg")

# Solve the problem
solver.solve()

# Retrieve the solution
v_opt = solver.getOutput("x")
u0_opt = [v_opt[actionIdx(0,k)] for k in range(nk)]
x0_opt = [v_opt[stateIdx(0,k)] for k in range(nk+1)]
x1_opt = [v_opt[stateIdx(1,k)] for k in range(nk+1)]
x2_opt = [v_opt[stateIdx(2,k)] for k in range(nk+1)]

# Retrieve the solution
v_opt = NP.array(solver.getOutput("x"))

# Get values at the beginning of each finite element
tgrid_u = NP.linspace(0,tf,nk)
tgrid_x = NP.linspace(0,tf,nk+1)


# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid_x,x0_opt,'--')
plt.plot(tgrid_x,x1_opt,'-')
plt.plot(tgrid_u,u0_opt,'-.')
plt.title("Van der Pol optimization - multiple shooting")
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

