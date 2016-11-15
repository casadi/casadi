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

# Declare variables (use simple, efficient DAG)
x0=SX.sym("x0"); x1=SX.sym("x1")
x = vertcat(x0,x1)

# Control
u = SX.sym("u")

# ODE right hand side
xdot = vertcat((1 - x1*x1)*x0 - x1 + u, x0)

# Lagrangian function
L = x0*x0 + x1*x1 + u*u

# Costate
lam = SX.sym("lam",2)

# Hamiltonian function
H = dot(lam,xdot) + L

# Costate equations
ldot = -gradient(H,x)

## The control must minimize the Hamiltonian, which is:
print("Hamiltonian: ", H)

# H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = gradient(H,u)     # this gives us 2*u + p
p = substitute(p,u,0) # replace u with zero: gives us p

# H's unconstrained minimizer is: u = -p/2
u_opt = -p/2

# We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = fmin(u_opt,1.0)
u_opt = fmax(u_opt,-0.75)
print("optimal control: ", u_opt)

# Augment f with lam_dot and subtitute in the value for the optimal control
f = vertcat(xdot,ldot)
f = substitute(f,u,u_opt)

# Function for obtaining the optimal control from the augmented state
u_fcn = Function("ufcn", [vertcat(x,lam)], [u_opt])

# Formulate the DAE
dae = {'x':vertcat(x,lam), 'ode':f}

# Augmented DAE state dimension
nX = 4

# End time
tf = 10.0

# Number of shooting nodes
num_nodes = 20

# Create an integrator (CVodes)
iopts = {}
iopts["abstol"] = 1e-8 # abs. tolerance
iopts["reltol"] = 1e-8 # rel. tolerance
iopts["t0"] = 0.0
iopts["tf"] = tf/num_nodes
I = integrator("I", "cvodes", dae, iopts)

# Variables in the root finding problem
NV = nX*(num_nodes+1)
V = MX.sym("V",NV)

# Get the state at each shooting node
X = []
v_offset = 0
for k in range(num_nodes+1):
  X.append(V[v_offset:v_offset+nX])
  v_offset = v_offset+nX

# Formulate the root finding problem
G = []
G.append(X[0][:2] - NP.array([0,1])) # states fixed, costates free at initial time
for k in range(num_nodes):
  XF = I(x0=X[k])["xf"]
  G.append(XF-X[k+1])
G.append(X[num_nodes][2:] - NP.array([0,0])) # costates fixed, states free at final time

# Terminal constraints: lam = 0
rfp = Function('rfp', [V], [vertcat(*G)])

# Select a solver for the root-finding problem
Solver = "nlpsol"
#Solver = "newton"
#Solver = "kinsol"

# Solver options
opts = {}
if Solver=="nlpsol":
    opts["nlpsol"] = "ipopt"
    opts["nlpsol_options"] = {"ipopt.hessian_approximation":"limited-memory"}
elif Solver=="newton":
    opts["linear_solver"] = CSparse
elif Solver=="kinsol":
    opts["linear_solver_type"] = "user_defined"
    opts["linear_solver"] = CSparse
    opts["max_iter"] = 1000

# Allocate a solver
solver = rootfinder('solver', Solver, rfp, opts)

# Solve the problem
V_sol = solver(0)

# Time grid for visualization
tgrid = NP.linspace(0,tf,100)

# Simulator to get optimal state and control trajectories
simulator = integrator('simulator', 'cvodes', dae, {'grid':tgrid,'output_t0':True})

# Simulate to get the trajectories
sol = simulator(x0 = V_sol[0:4])["xf"]

# Calculate the optimal control
ufcn_all = u_fcn.map(len(tgrid))
u_opt = ufcn_all(sol)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid, sol[0, :].T, '--')
plt.plot(tgrid, sol[1, :].T, '-')
plt.plot(tgrid, u_opt.T, '-.')
plt.title("Van der Pol optimization - indirect multiple shooting")
plt.xlabel('time')
plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.grid()
plt.show()
