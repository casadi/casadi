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

# For documentation about this examples, check http://docs.casadi.org/documents/mhe_spring_damper.pdf

from casadi import *
import numpy as NP
import matplotlib.pyplot as plt
import time
from casadi.tools import *
from scipy import linalg, matrix
plt.interactive(True)

NP.random.seed(0)

# Settings of the filter
N = 10 # Horizon length
dt = 0.05; # Time step

sigma_p = 0.005 # Standard deviation of the position measurements
sigma_w = 0.1 # Standard deviation for the process noise
R = DM(1/sigma_p**2) # resulting weighting matrix for the position measurements
Q = DM(1/sigma_w**2) # resulting weighting matrix for the process noise

Nsimulation = 1000 # Lenght of the simulation

# Parameters of the system
m = 1 # The weight of the mass
k = 1 # The spring constant
c = 0.5 # The damping of the system
# The state
states = struct_symSX(["x","dx"]) # Full state vector of the system: position x and velocity dx
Nstates = states.size # Number of states
# Set up some aliases
x,dx = states[...]

# The control input
controls = struct_symSX(["F"]) # Full control vector of the system: Input force F
Ncontrols = controls.size # Number of control inputs
# Set up some aliases
F, = controls[...]

# Disturbances
disturbances = struct_symSX(["w"]) # Process noise vector
Ndisturbances = disturbances.size # Number of disturbances
# Set up some aliases
w, = disturbances[...]

# Measurements
measurements = struct_symSX(["y"]) # Measurement vector
Nmeas = measurements.size # Number of measurements
# Set up some aliases
y, = measurements[...]

# Create Structure for the entire horizon

# Structure that will be degrees of freedom for the optimizer
shooting = struct_symSX([(entry("X",repeat=N,struct=states),entry("W",repeat=N-1,struct=disturbances))])
# Structure that will be fixed parameters for the optimizer
parameters = struct_symSX([(entry("U",repeat=N-1,struct=controls),entry("Y",repeat=N,struct=measurements),entry("S",shape=(Nstates,Nstates)),entry("x0",shape=(Nstates,1)))])
S = parameters["S"]
x0 = parameters["x0"]
# Define the ODE right hand side
rhs = struct_SX(states)
rhs["x"] = dx
rhs["dx"] = (-k*x-c*dx+F)/m+w

f = Function('f', [states,controls,disturbances],[rhs])

# Build an integrator for this system: Runge Kutta 4 integrator
k1 = f(states,controls,disturbances)
k2 = f(states+dt/2.0*k1,controls,disturbances)
k3 = f(states+dt/2.0*k2,controls,disturbances)
k4 = f(states+dt*k3,controls,disturbances)

states_1 = states+dt/6.0*(k1+2*k2+2*k3+k4)
phi = Function('phi', [states,controls,disturbances],[states_1])
PHI = phi.jacobian_old(0, 0)
# Define the measurement system
h = Function('h', [states],[x]) # We have measurements of the position
H = h.jacobian_old(0, 0)
# Build the objective
obj = 0
# First the arrival cost
obj += mtimes([(shooting["X",0]-parameters["x0"]).T,S,(shooting["X",0]-parameters["x0"])])
#Next the cost for the measurement noise
for i in range(N):
  vm = h(shooting["X",i])-parameters["Y",i]
  obj += mtimes([vm.T,R,vm])
#And also the cost for the process noise
for i in range(N-1):
  obj += mtimes([shooting["W",i].T,Q,shooting["W",i]])

# Build the multiple shooting constraints
g = []
for i in range(N-1):
  g.append( shooting["X",i+1] - phi(shooting["X",i],parameters["U",i],shooting["W",i]) )

# Formulate the NLP
nlp = {'x':shooting, 'p':parameters, 'f':obj, 'g':vertcat(*g)}

# Make a simulation to create the data for the problem
simulated_X = DM.zeros(Nstates,Nsimulation)
simulated_X[:,0] = DM([1,0]) # Initial state
t = NP.linspace(0,(Nsimulation-1)*dt,Nsimulation) # Time grid
simulated_U = DM(cos(t[0:-1])).T # control input for the simulation
simulated_U[:,int(Nsimulation/2):] = 0.0
simulated_W = DM(sigma_w*NP.random.randn(Ndisturbances,Nsimulation-1)) # Process noise for the simulation
for i in range(Nsimulation-1):
  simulated_X[:,i+1] = phi(simulated_X[:,i], simulated_U[:,i], simulated_W[:,i])
#Create the measurements from these states
simulated_Y = DM.zeros(Nmeas,Nsimulation) # Holder for the measurements
for i in range(Nsimulation):
  simulated_Y[:,i] = h(simulated_X[:,i])
# Add noise the the position measurements
simulated_Y += sigma_p*NP.random.randn(simulated_Y.shape[0],simulated_Y.shape[1])

#The initial estimate and related covariance, which will be used for the arrival cost
sigma_x0 = 0.01
P = sigma_x0**2*DM.eye(Nstates)
x0 = simulated_X[:,0] + sigma_x0*NP.random.randn(Nstates,1)
# Create the solver
opts = {"ipopt.print_level":0, "print_time": False, 'ipopt.max_iter':100}
nlpsol = nlpsol("nlpsol", "ipopt", nlp, opts)

# Create a holder for the estimated states and disturbances
estimated_X= DM.zeros(Nstates,Nsimulation)
estimated_W = DM.zeros(Ndisturbances,Nsimulation-1)
# For the first instance we run the filter, we need to initialize it.
current_parameters = parameters(0)
current_parameters["U",lambda x: horzcat(*x)] = simulated_U[:,0:N-1]
current_parameters["Y",lambda x: horzcat(*x)] = simulated_Y[:,0:N]
current_parameters["S"] = linalg.inv(P) # Arrival cost is the inverse of the initial covariance
current_parameters["x0"] = x0
initialisation_state = shooting(0)
initialisation_state["X",lambda x: horzcat(*x)] = simulated_X[:,0:N]
res = nlpsol(p=current_parameters, x0=initialisation_state, lbg=0, ubg=0)

# Get the solution
solution = shooting(res["x"])
estimated_X[:,0:N] = solution["X",lambda x: horzcat(*x)]
estimated_W[:,0:N-1] = solution["W",lambda x: horzcat(*x)]

# Now make a loop for the rest of the simulation
for i in range(1,Nsimulation-N+1):

  # Update the arrival cost, using linearisations around the estimate of MHE at the beginning of the horizon (according to the 'Smoothed EKF Update'): first update the state and covariance with the measurement that will be deleted, and next propagate the state and covariance because of the shifting of the horizon
  print("step %d/%d (%s)" % (i, Nsimulation-N , nlpsol.stats()["return_status"]))
  H0 = H(solution["X",0])[0]
  K = mtimes([P,H0.T,linalg.inv(mtimes([H0,P,H0.T])+R)])
  P = mtimes((DM.eye(Nstates)-mtimes(K,H0)),P)
  h0 = h(solution["X",0])
  x0 = x0 + mtimes(K, current_parameters["Y",0]-h0-mtimes(H0,x0-solution["X",0]))
  x0 = phi(x0, current_parameters["U",0], solution["W",0])
  F = PHI(solution["X",0], current_parameters["U",0], solution["W",0])[0]
  P = mtimes([F,P,F.T]) + linalg.inv(Q)
  # Get the measurements and control inputs
  current_parameters["U",lambda x: horzcat(*x)] = simulated_U[:,i:i+N-1]
  current_parameters["Y",lambda x: horzcat(*x)] = simulated_Y[:,i:i+N]
  current_parameters["S"] = linalg.inv(P)
  current_parameters["x0"] = x0
  # Initialize the system with the shifted solution
  initialisation_state["W",lambda x: horzcat(*x),0:N-2] = estimated_W[:,i:i+N-2] # The shifted solution for the disturbances
  initialisation_state["W",N-2] = DM.zeros(Ndisturbances,1) # The last node for the disturbances is initialized with zeros
  initialisation_state["X",lambda x: horzcat(*x),0:N-1] = estimated_X[:,i:i+N-1] # The shifted solution for the state estimates
  # The last node for the state is initialized with a forward simulation
  phi0 = phi(initialisation_state["X",N-1], current_parameters["U",-1], initialisation_state["W",-1])
  initialisation_state["X",N-1] = phi0
  # And now initialize the solver and solve the problem
  res = nlpsol(p=current_parameters, x0=initialisation_state, lbg=0, ubg=0)
  solution = shooting(res["x"])

  # Now get the state estimate. Note that we are only interested in the last node of the horizon
  estimated_X[:,N-1+i] = solution["X",N-1]
  estimated_W[:,N-2+i] = solution["W",N-2]
# Plot the results
plt.figure(1)
plt.clf()
plt.plot(t,vec(estimated_X[0,:]),'b--')
plt.plot(t,vec(simulated_X[0,:]),'r--')
plt.title("Position")
plt.xlabel('Time')
plt.legend(['Estimated position','Real position'])
plt.grid()

plt.figure(2)
plt.clf()
plt.plot(t,vec(estimated_X[0,:]-simulated_X[0,:]),'b--')
plt.title("Position error")
plt.xlabel('Time')
plt.legend(['Error between estimated and real position'])
plt.grid()

plt.show()

error = estimated_X[0,:]-simulated_X[0,:]
print(mtimes(error,error.T))
assert(mtimes(error,error.T)<0.01)
