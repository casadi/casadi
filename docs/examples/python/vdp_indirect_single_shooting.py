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
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt


# Declare variables (use simple, efficient DAG)
x0=SX.sym('x0'); x1=SX.sym('x1')
x = vertcat(x0,x1)

# Control
u = SX.sym('u')

# ODE right hand side
xdot = vertcat((1 - x1*x1)*x0 - x1 + u, x0)

# Lagrangian function
L = x0*x0 + x1*x1 + u*u

# Costate
lam = SX.sym('lam',2)

# Hamiltonian function
H = dot(lam,xdot) + L

# Costate equations
ldot = -gradient(H,x)

## The control must minimize the Hamiltonian, which is:
print('Hamiltonian: ', H)

# H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = gradient(H,u)    # this gives us 2*u + p
p = substitute(p,u,0) # replace u with zero: gives us p

# H's unconstrained minimizer is: u = -p/2
u_opt = -p/2

# We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = fmin(u_opt, 1.0)
u_opt = fmax(u_opt, -0.75)
print('optimal control: ', u_opt)

# Augment f with lam_dot and substitute in the value for the optimal control
f = vertcat(xdot,ldot)
f = substitute(f,u,u_opt)

# Function for obtaining the optimal control from the augmented state
u_fcn = Function('ufcn', [vertcat(x,lam)], [u_opt])

# Formulate the DAE
dae = {'x':vertcat(x,lam), 'ode':f}

# Create an integrator (CVodes)
opts = {}
opts['abstol'] = 1e-8 # abs. tolerance
opts['reltol'] = 1e-8 # rel. tolerance
I = integrator('I', 'cvodes', dae, 0, 10, opts)

# The initial state
x_init = NP.array([0.,1.])

# The initial costate
l_init = MX.sym('l_init',2)

# The initial condition for the shooting
X = vertcat(x_init, l_init)

# Call the integrator
X = I(x0=X)['xf']

# Costate at the final time should be zero (cf. Bryson and Ho)
lam_f = X[2:4]
g = lam_f

# Formulate root-finding problem
rfp = Function('rfp', [l_init], [g])

# Select a solver for the root-finding problem
Solver = 'nlpsol'
#Solver = 'newton'
#Solver = 'kinsol'

# Allocate an implict solver
opts = {}
if Solver=='nlpsol':
    opts['nlpsol'] = 'ipopt'
    opts['nlpsol_options'] = {'ipopt.hessian_approximation':'limited-memory'}
elif Solver=='kinsol':
    opts['linear_solver_type'] = 'user_defined'
    opts['max_iter'] = 1000
solver = rootfinder('solver', Solver, rfp, opts)

# Solve the problem
l_init_opt = solver(0)
l_init_opt = NP.array(l_init_opt.nonzeros())

# Time grid for visualization
tgrid = NP.linspace(0, 10, 100)

# Simulator to get optimal state and control trajectories
simulator = integrator('simulator', 'cvodes', dae, 0, tgrid)

# Simulate to get the state trajectory
sol = simulator(x0 = NP.concatenate((x_init, l_init_opt)))['xf']

# Calculate the optimal control
ufcn_all = u_fcn.map(len(tgrid))
u_opt = ufcn_all(sol)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid, sol[0, :].T, '--')
plt.plot(tgrid, sol[1, :].T, '-')
plt.plot(tgrid, u_opt.T, '-.')
plt.title('Van der Pol optimization - indirect single shooting')
plt.xlabel('time')
plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.grid()
plt.show()
