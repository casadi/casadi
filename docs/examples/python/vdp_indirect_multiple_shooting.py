#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
import casadi as ca
import numpy as NP
import matplotlib.pyplot as plt

# Declare variables (use simple, efficient DAG)
x0=ca.SX.sym("x0"); x1=ca.SX.sym("x1")
x = ca.vertcat(x0,x1)

# Control
u = ca.SX.sym("u")

# ODE right hand side
xdot = ca.vertcat((1 - x1*x1)*x0 - x1 + u, x0)

# Lagrangian function
L = x0*x0 + x1*x1 + u*u

# Costate
lam = ca.SX.sym("lam",2)

# Hamiltonian function
H = ca.dot(lam,xdot) + L

# Costate equations
ldot = -ca.gradient(H,x)

## The control must minimize the Hamiltonian, which is:
print("Hamiltonian: ", H)

# H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = ca.gradient(H,u)     # this gives us 2*u + p
p = ca.substitute(p,u,0) # replace u with zero: gives us p

# H's unconstrained minimizer is: u = -p/2
u_opt = -p/2

# We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = ca.fmin(u_opt,1.0)
u_opt = ca.fmax(u_opt,-0.75)
print("optimal control: ", u_opt)

# Augment f with lam_dot and subtitute in the value for the optimal control
f = ca.vertcat(xdot,ldot)
f = ca.substitute(f,u,u_opt)

# Function for obtaining the optimal control from the augmented state
u_fcn = ca.Function("ufcn", [ca.vertcat(x,lam)], [u_opt])

# Formulate the DAE
dae = {'x':ca.vertcat(x,lam), 'ode':f}

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
I = ca.integrator("I", "cvodes", dae, 0, tf/num_nodes, iopts)

# Variables for the states at each shooting node
X = ca.MX.sym('X',nX,num_nodes+1)

# Formulate the root finding problem
G = []
G.append(X[:2,0] - ca.vertcat(0,1)) # states fixed, costates free at initial time
for k in range(num_nodes):
  XF = I(x0=X[:,k])["xf"]
  G.append(XF-X[:,k+1])
G.append(X[2:,num_nodes] - ca.vertcat(0,0)) # costates fixed, states free at final time

# Terminal constraints: lam = 0
rfp = {"x": ca.vec(X), "g": ca.vertcat(*G)}

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
    opts["linear_solver"] = "csparse"
elif Solver=="kinsol":
    opts["linear_solver_type"] = "user_defined"
    opts["linear_solver"] = "csparse"
    opts["max_iter"] = 1000

# Allocate a solver
solver = ca.rootfinder('solver', Solver, rfp, opts)

# Solve the problem
X_sol = solver(x0=0)['x']
print(X_sol)

# Time grid for visualization
tgrid = NP.linspace(0,tf,100)

# Simulator to get optimal state and control trajectories
simulator = ca.integrator('simulator', 'cvodes', dae, 0, tgrid)

# Simulate to get the trajectories
sol = simulator(x0 = X_sol[0:4])["xf"]

# Calculate the optimal control
u_opt = u_fcn(sol)

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

print(sol[0,:])
