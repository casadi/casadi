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

# For documentation about this examples, check http://docs.casadi.org/documents/mhe_spring_damper.pdf

from casadi import *
import numpy as NP
import matplotlib.pyplot as plt
import time
from casadi.tools import *
from scipy import linalg, matrix
plt.interactive(True)

# Settings of the filter
N = 10 # Horizon length
dt = 0.05; # Time step

sigma_p = 0.005 # Standard deviation of the position measurements
sigma_w = 0.1 # Standard deviation for the process noise
R = DMatrix(1/sigma_p**2) # resulting weighting matrix for the position measurements
Q = DMatrix(1/sigma_w**2) # resulting weighting matrix for the process noise

Nsimulation = 1000 # Lenght of the simulation

# Parameters of the system
m = 1 # The weight of the mass
k = 1 # The spring constant
c = 0.5 # The damping of the system
# The state
states = struct_ssym(["x","dx"]) # Full state vector of the system: position x and velocity dx
Nstates = states.size # Number of states
# Set up some aliases
x,dx = states[...]

# The control input
controls = struct_ssym(["F"]) # Full control vector of the system: Input force F
Ncontrols = controls.size # Number of control inputs
# Set up some aliases
F, = controls[...]

# Disturbances
disturbances = struct_ssym(["w"]) # Process noise vector
Ndisturbances = disturbances.size # Number of disturbances
# Set up some aliases
w, = disturbances[...]

# Measurements
measurements = struct_ssym(["y"]) # Measurement vector
Nmeas = measurements.size # Number of measurements
# Set up some aliases
y, = measurements[...]

# Create Structure for the entire horizon

# Structure that will be degrees of freedom for the optimizer
shooting = struct_ssym([(entry("X",repeat=N,struct=states),entry("W",repeat=N-1,struct=disturbances))])
# Structure that will be fixed parameters for the optimizer
parameters = struct_ssym([(entry("U",repeat=N-1,struct=controls),entry("Y",repeat=N,struct=measurements),entry("S",shape=(Nstates,Nstates)),entry("x0",shape=(Nstates,1)))])
S = parameters["S"]
x0 = parameters["x0"]
# Define the ODE right hand side
rhs = struct_SX(states)
rhs["x"] = dx
rhs["dx"] = (-k*x-c*dx+F)/m+w

f = SXFunction([states,controls,disturbances],[rhs])
f.init()

# Build an integrator for this system: Runge Kutta 4 integrator
k1 = f.eval([states,controls,disturbances])[0]
k2 = f.eval([states+dt/2.0*k1,controls,disturbances])[0]
k3 = f.eval([states+dt/2.0*k2,controls,disturbances])[0]
k4 = f.eval([states+dt*k3,controls,disturbances])[0]

states_1 = states+dt/6.0*(k1+2*k2+2*k3+k4)
phi = SXFunction([states,controls,disturbances],[states_1])
phi.init()
PHI = phi.jacobian()
PHI.init()
# Define the measurement system
h = SXFunction([states],[x]) # We have measurements of the position
h.init()
H = h.jacobian()
H.init()
# Build the objective
obj = 0
# First the arrival cost
obj += mul([(shooting["X",0]-parameters["x0"]).T,S,(shooting["X",0]-parameters["x0"])])
#Next the cost for the measurement noise
for i in range(N):
  vm = h.eval([shooting["X",i]])[0]-parameters["Y",i]
  obj += mul([vm.T,R,vm])
#And also the cost for the process noise
for i in range(N-1):
  obj += mul([shooting["W",i].T,Q,shooting["W",i]])

f_obj = SXFunction([shooting,parameters],[obj])
f_obj.init()
# Build the multiple shooting constraints
g = []
for i in range(N-1):
  g.append( shooting["X",i+1] - phi.eval([shooting["X",i],parameters["U",i],shooting["W",i]])[0] )
G = SXFunction([shooting,parameters],[vertcat(g)])
G.init()
# Make a simulation to create the data for the problem
simulated_X = DMatrix.zeros(Nstates,Nsimulation)
simulated_X[:,0] = DMatrix([1,0]) # Initial state
t = NP.linspace(0,(Nsimulation-1)*dt,Nsimulation) # Time grid
simulated_U = DMatrix(cos(t[0:-1])).T # control input for the simulation
simulated_U[:,Nsimulation/2:] = 0.0
simulated_W = DMatrix(sigma_w*NP.random.randn(Ndisturbances,Nsimulation-1)) # Process noise for the simulation
for i in range(Nsimulation-1):
  phi.input(0).set(simulated_X[:,i])
  phi.input(1).set(simulated_U[:,i])
  phi.input(2).set(simulated_W[:,i])
  phi.evaluate()
  simulated_X[:,i+1] = phi.output(0)
#Create the measurements from these states
simulated_Y = DMatrix(Nmeas,Nsimulation) # Holder for the measurements
for i in range(Nsimulation):
  h.input(0).set(simulated_X[:,i])
  h.evaluate()
  simulated_Y[:,i] = h.output(0)
# Add noise the the position measurements
simulated_Y += sigma_p*NP.random.randn(simulated_Y.shape[0],simulated_Y.shape[1])

#The initial estimate and related covariance, which will be used for the arrival cost
sigma_x0 = 0.01
P = sigma_x0**2*DMatrix.eye(Nstates)
x0 = simulated_X[:,0] + sigma_x0*NP.random.randn(Nstates,1)
# Create the solver
nlp_solver = IpoptSolver(f_obj,G)
nlp_solver.setOption({"print_level":0, "print_time": False})
nlp_solver.setOption('linear_solver','MA57')
nlp_solver.setOption('max_iter',100)
nlp_solver.init()

# Set the bounds for the constraints: we only have the multiple shooting constraints, so all constraints have upper and lower bound of zero
nlp_solver.input("lbg").setAll(0)
nlp_solver.input("ubg").setAll(0)

# Create a holder for the estimated states and disturbances
estimated_X= DMatrix(Nstates,Nsimulation)
estimated_W = DMatrix(Ndisturbances,Nsimulation-1)
# For the first instance we run the filter, we need to initialize it.
current_parameters = parameters(0)
current_parameters["U",horzcat] = simulated_U[:,0:N-1]
current_parameters["Y",horzcat] = simulated_Y[:,0:N]
current_parameters["S"] = linalg.inv(P) # Arrival cost is the inverse of the initial covariance
current_parameters["x0"] = x0
initialisation_state = shooting(0)
initialisation_state["X",horzcat] = simulated_X[:,0:N]

nlp_solver.setInput(current_parameters,"p")
nlp_solver.setInput(initialisation_state,"x0")

nlp_solver.solve()
# Get the solution
solution = shooting(nlp_solver.output("x"))
estimated_X[:,0:N] = solution["X",horzcat]
estimated_W[:,0:N-1] = solution["W",horzcat]

# Now make a loop for the rest of the simulation
for i in range(1,Nsimulation-N+1):
  # Update the arrival cost, using linearisations around the estimate of MHE at the beginning of the horizon (according to the 'Smoothed EKF Update'): first update the state and covariance with the measurement that will be deleted, and next propagate the state and covariance because of the shifting of the horizon
  print "step %d/%d (%s)" % (i, Nsimulation-N , nlp_solver.getStat("return_status"))
  H.input(0).set(solution["X",0])
  H.evaluate()
  H0 = H.output(0)
  K = mul([P,H0.T,linalg.inv(mul([H0,P,H0.T])+R)])
  P = mul((DMatrix.eye(Nstates)-mul(K,H0)),P)
  h.input(0).set(solution["X",0])
  h.evaluate()
  x0 = x0 + mul(K, current_parameters["Y",0]-h.output(0)-mul(H0,x0-solution["X",0]) )
  phi.input(0).set(x0)
  phi.input(1).set(current_parameters["U",0])
  phi.input(2).set(solution["W",0])
  phi.evaluate()
  x0 = phi.output(0)
  PHI.input(0).set(solution["X",0])
  PHI.input(1).set(current_parameters["U",0])
  PHI.input(2).set(solution["W",0])
  PHI.evaluate()
  F = PHI.output(0)
  PHI.evaluate()
  P = mul([F,P,F.T]) + linalg.inv(Q)
  # Get the measurements and control inputs 
  current_parameters["U",horzcat] = simulated_U[:,i:i+N-1]
  current_parameters["Y",horzcat] = simulated_Y[:,i:i+N]
  current_parameters["S"] = linalg.inv(P)
  current_parameters["x0"] = x0
  # Initialize the system with the shifted solution
  initialisation_state["W",horzcat,0:N-2] = estimated_W[:,i:i+N-2] # The shifted solution for the disturbances
  initialisation_state["W",N-2] = DMatrix.zeros(Ndisturbances,1) # The last node for the disturbances is initialized with zeros
  initialisation_state["X",horzcat,0:N-1] = estimated_X[:,i:i+N-1] # The shifted solution for the state estimates
  # The last node for the state is initialized with a forward simulation
  phi.input(0).set(initialisation_state["X",N-1] )
  phi.input(1).set(current_parameters["U",-1])
  phi.input(2).set(initialisation_state["W",-1])
  phi.evaluate()
  initialisation_state["X",N-1] = phi.output(0)
  # And now initialize the solver and solve the problem
  nlp_solver.setInput(current_parameters,"p")
  nlp_solver.setInput(initialisation_state,"x0")
  nlp_solver.solve()
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
print mul(error,error.T)
