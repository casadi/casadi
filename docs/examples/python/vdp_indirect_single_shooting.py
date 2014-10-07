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
x = vertcat((x0,x1))

# Control
u = SX.sym("u")

# ODE right hand side
xdot = vertcat([(1 - x1*x1)*x0 - x1 + u, x0])

# Lagrangian function
L = x0*x0 + x1*x1 + u*u

# Costate
lam = SX.sym("lam",2)

# Hamiltonian function
H = inner_prod(lam,xdot) + L

# Costate equations
ldot = -gradient(H,x)

## The control must minimize the Hamiltonian, which is:
print "Hamiltonian: ", H

# H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = gradient(H,u)    # this gives us 2*u + p
p = substitute(p,u,0) # replace u with zero: gives us p

# H's unconstrained minimizer is: u = -p/2
u_opt = -p/2

# We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = min(u_opt,1.0)
u_opt = max(u_opt,-0.75)
print "optimal control: ", u_opt

# Augment f with lam_dot and substitute in the value for the optimal control
f = vertcat((xdot,ldot))
f = substitute(f,u,u_opt)

# Create the right hand side function
rhs_in = daeIn(x=vertcat((x,lam)))
rhs = SXFunction(rhs_in,daeOut(ode=f))

# Create an integrator (CVodes)
I = Integrator("cvodes", rhs)
I.setOption("abstol",1e-8) # abs. tolerance
I.setOption("reltol",1e-8) # rel. tolerance
I.setOption("t0",0.0)
I.setOption("tf",10.0)
I.init()

# The initial state
x_init = NP.array([0.,1.])

# The initial costate
l_init = MX.sym("l_init",2)

# The initial condition for the shooting
X = vertcat((x_init,l_init))

# Call the integrator
X, = integratorOut(I.call(integratorIn(x0=X)),"xf")

# Costate at the final time should be zero (cf. Bryson and Ho)
lam_f = X[2:4]
g = lam_f

# Formulate root-finding problem
rfp = MXFunction([l_init],[g])

# Select a solver for the root-finding problem
Solver = "nlp"
#Solver = "newton"
#Solver = "kinsol"

# Allocate an implict solver
solver = ImplicitFunction(Solver, rfp)
if Solver=="nlp":
    solver.setOption("nlp_solver", "ipopt")
    solver.setOption("nlp_solver_options",{"hessian_approximation":"limited-memory"})
elif Solver=="newton":
    solver.setOption("linear_solver",CSparse)
elif Solver=="kinsol":
    solver.setOption("linear_solver_type","user_defined")
    solver.setOption("linear_solver",CSparse)
    solver.setOption("max_iter",1000)

# Initialize the solver
solver.init()

# Pass initial guess
#solver.setInput([   0,   0], "x0")

# Solve the problem
solver.evaluate()

# Retrieve the optimal solution
l_init_opt = NP.array(solver.output().data())

# Time grid for visualization
tgrid = NP.linspace(0,10,100)

# Output functions
output_fcn = SXFunction(rhs_in,[x0,x1,u_opt])

# Simulator to get optimal state and control trajectories
simulator = Simulator(I, output_fcn, tgrid)
simulator.init()

# Pass initial conditions to the simulator
simulator.setInput(NP.concatenate((x_init,l_init_opt)),"x0")

# Simulate to get the trajectories
simulator.evaluate()

# Get optimal control
x_opt = simulator.getOutput(0).T
y_opt = simulator.getOutput(1).T
u_opt = simulator.getOutput(2).T

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid,x_opt,'--')
plt.plot(tgrid,y_opt,'-')
plt.plot(tgrid,u_opt,'-.')
plt.title("Van der Pol optimization - indirect single shooting")
plt.xlabel('time')
plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.grid()
plt.show()
