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

nk = 60      # Control discretization
coll = True # Use collocation integrator
coll = False

# Declare variables (use scalar graph)
numStates = 5
numActions = 1
numParams = 0

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

def paramIdx(paramK):
  if paramK > numParams - 1:
    raise ValueError("paramK > numParams - 1")
  if paramK < 0:
    raise ValueError("paramK < 0")
  return numStates*(nk+1) + numActions*nk + paramK

# Total number of variables
nv = numActions*nk + numStates*(nk+1) + numParams

# Declare variable vector
V = msym("V", nv)

# end time
#tf = V[paramIdx(0)]
tf = 8.0

t  = ssym("t")    # time
u  = ssym("u",numActions) # control
x  = ssym("x",numStates)  # state
xp = ssym("xd",numStates) # state derivative

# ODE/DAE residual function
mass = 1.0
gravity = 9.8
def basketDxdt():
    q  = x[0]
    qd = x[1]
    L  = x[2]
    Ld = x[3]

    T = u[0]

    g = gravity
    m = mass

    qdd = (-g*sin(q) - 2*Ld*qd)/L
    Ldd = g*cos(q) + L*qd*qd - T/m

#    power = -m*g*sin(q)*L*qd + (m*g*cos(q)-T)*Ld
    power = -T*Ld

    return [qd,
            qdd,
            Ld,
            Ldd,
            power]

res = basketDxdt() - xp
f = SXFunction([t,x,u,xp],[res])
f.init()

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
  f_d.setOption("abstol",1e-9) # tolerance
  f_d.setOption("reltol",1e-6) # tolerance
  f_d.setOption("steps_per_checkpoint",1000)

f_d.setOption("tf",tf/nk) # final time
f_d.init()

# Variable bounds initialized to +/- inf
VMIN = -inf*NP.ones(nv)
VMAX = inf*NP.ones(nv)

## Control bounds
for k in range(nk):
  # tension
  VMIN[actionIdx(0,k)] = 0.0
  VMAX[actionIdx(0,k)] = 20.0

# State bounds
for k in range(nk+1):
  # theta
  VMIN[stateIdx(0,k)] = -0.5*pi
  VMAX[stateIdx(0,k)] = 0.5*pi

  # thetaDot
  VMIN[stateIdx(1,k)] = -5.0*pi
  VMAX[stateIdx(1,k)] = 5.0*pi

  # length
  VMIN[stateIdx(2,k)] = 0.5
  VMAX[stateIdx(2,k)] = 30.0

  # lengthDot
  VMIN[stateIdx(3,k)] = -100.0
  VMAX[stateIdx(3,k)] = 100.0

# param bounds
if numParams > 0:
  VMIN[paramIdx(0)] = 9.0
  VMAX[paramIdx(0)] = 11.0


# Initial condition
VMIN[stateIdx(0,0)] = VMAX[stateIdx(0,0)] = 0.8
VMIN[stateIdx(1,0)] = VMAX[stateIdx(1,0)] = 0
VMIN[stateIdx(2,0)] = VMAX[stateIdx(2,0)] = 10
VMIN[stateIdx(3,0)] = VMAX[stateIdx(3,0)] = 0
VMIN[stateIdx(4,0)] = VMAX[stateIdx(4,0)] = 0

# Terminal constraint
VMIN[stateIdx(3,nk)] = VMAX[stateIdx(3,nk)] = 0
# theta
VMIN[stateIdx(0,nk)] = -0.3
VMAX[stateIdx(0,nk)] =  0.3
# thetaDot
VMIN[stateIdx(1,nk)] = -0.4
VMAX[stateIdx(1,nk)] =  0.4

# length
VMAX[stateIdx(2,nk)] =  1


# Initial solution guess
VINIT = NP.zeros(nv)
# initial state
for k in range(nk+1):
  # theta
  VINIT[stateIdx(0,k)] = 0
  # thetaDot
  VINIT[stateIdx(1,k)] = 0
  # length
  VINIT[stateIdx(2,k)] = 10
  # lengthDot
  VINIT[stateIdx(3,k)] = 0
# initial action
for k in range(nk):
  VINIT[actionIdx(0,k)] = 0
# initial params
if numParams > 0:
  VINIT[paramIdx(0)] = 10.0

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
def totalEnergy(vec, k):
  q  = vec[stateIdx(0,k)]
  qd = vec[stateIdx(1,k)]
  L  = vec[stateIdx(2,k)]
  Ld = vec[stateIdx(3,k)]
  ke = 0.5*mass*(L*L*qd*qd + Ld*Ld)
  pe = -mass*gravity*L*cos(q)
  return ke + pe
F = MXFunction([V],[totalEnergy(V,nk)])

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
v_opt = NP.array(solver.getOutput("x"))

# Get values at the beginning of each finite element
tgrid_u = NP.linspace(0,tf,nk)
tgrid_x = NP.linspace(0,tf,nk+1)

# Plot the results
def makePlots(blah):
  plt.figure(1)
  plt.clf()
  for k,boo in enumerate(blah):
    plt.subplot(len(blah), 1, k+1)
    plt.plot(boo[0], boo[1])
    plt.ylabel(boo[2])

  plt.xlabel('time')
  plt.grid()
  plt.show()

makePlots([(tgrid_x, [v_opt[stateIdx(0,k)] for k in range(nk+1)], 'q'),
           (tgrid_x, [v_opt[stateIdx(1,k)] for k in range(nk+1)], 'dq/dt'),
           (tgrid_x, [v_opt[stateIdx(2,k)] for k in range(nk+1)], 'L'),
           (tgrid_x, [v_opt[stateIdx(3,k)] for k in range(nk+1)], 'dL/dt'),
           (tgrid_x, [v_opt[stateIdx(4,k)] for k in range(nk+1)], 'int(power)'),
           (tgrid_x, [totalEnergy(v_opt,k)-totalEnergy(v_opt,0) for k in range(nk+1)], 'KE+PE'),
           (tgrid_u, [-v_opt[stateIdx(3,k)]*v_opt[actionIdx(0,k)] for k in range(nk)], 'tension power'),
           (tgrid_u, [v_opt[actionIdx(0,k)] for k in range(nk)], 'Tension')])
