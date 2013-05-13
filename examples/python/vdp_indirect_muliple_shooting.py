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
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

# time
t = ssym("t")

# Declare variables (use simple, efficient DAG)
x0=ssym("x0"); x1=ssym("x1")
x = vertcat((x0,x1))

# Control
u = ssym("u")

# ODE right hand side
xdot = vertcat([(1 - x1*x1)*x0 - x1 + u, x0])

# Lagrangian function
L = x0*x0 + x1*x1 + u*u

# Costate
lam = ssym("lam",2)

# Hamiltonian function
H = inner_prod(lam,xdot) + L
Hfcn = SXFunction([x,lam,u,t],[H])
Hfcn.init()

# Costate equations
ldot = -Hfcn.grad(0,0)

## The control must minimize the Hamiltonian, which is:
print "Hamiltonian: ", H

# H is of a convect quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = Hfcn.grad(2,0)    # this gives us 2*u + p
p = substitute(p,u,0) # replace u with zero: gives us p

# H's unconstrained minimizer is: u = -p/2
u_opt = -p/2

# We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = min(u_opt,1.0)
u_opt = max(u_opt,-0.75)
print "optimal control: ", u_opt

# Augment f with lam_dot and subtitute in the value for the optimal control
f = vertcat((xdot,ldot))
f = substitute(f,u,u_opt)

# Create the right hand side function
rhs_in = list(daeIn(x=vertcat((x,lam))))
rhs_in[DAE_T] = t
rhs = SXFunction(rhs_in,daeOut(ode=f))

# Augmented DAE state dimension
nX = 4

# End time
tf = 10.0

# Number of shooting nodes
num_nodes = 20

# Create an integrator (CVodes)
I = CVodesIntegrator(rhs)
I.setOption("abstol",1e-8) # abs. tolerance
I.setOption("reltol",1e-8) # rel. tolerance
I.setOption("t0",0.0)
I.setOption("tf",tf/num_nodes)
I.init()

# Variables in the root finding problem
NV = nX*(num_nodes+1)
V = msym("V",NV)

# Get the state at each shooting node
X = []
v_offset = 0
for k in range(num_nodes+1):
  X.append(V[v_offset:v_offset+nX])
  v_offset = v_offset+nX

# Formulate the root finding problem
G = []
for k in range(num_nodes):
  XF, = integratorOut(I.call(integratorIn(x0=X[k])),"xf")
  G.append(XF-X[k+1])

# Terminal constraints: lam = 0
G = MXFunction([V],[vertcat(G)])

# Dummy objective function (there are no degrees of freedom)
F = MXFunction([V],[0])

# Allocate NLP solver
solver = IpoptSolver(F,G)

# Initialize the NLP solver
solver.init()

# Set bounds and initial guess
solver.setInput([   0,   1,-inf,-inf] + (num_nodes-1)*[-inf,-inf,-inf,-inf] + [-inf,-inf,   0,   0], "lbx")
solver.setInput([   0,   1, inf, inf] + (num_nodes-1)*[ inf, inf, inf, inf] + [ inf, inf,   0,   0], "ubx")
solver.setInput([   0,   0,   0,   0] + (num_nodes-1)*[   0,   0,   0,   0] + [   0,   0,   0,   0], "x0")
solver.setInput(num_nodes*[0,0,0,0], "lbg")
solver.setInput(num_nodes*[0,0,0,0], "ubg")

# Solve the problem
solver.solve()

# Time grid for visualization
tgrid = NP.linspace(0,tf,100)

# Output functions
output_fcn = SXFunction(rhs_in,[x0,x1,u_opt])

# Increase the end time for the integrator
I.setOption("tf",tf)
I.init()

# Simulator to get optimal state and control trajectories
simulator = Simulator(I, output_fcn, tgrid)
simulator.init()

# Pass initial conditions to the simulator
simulator.setInput(solver.output("x")[0:4],"x0")

# Simulate to get the trajectories
simulator.evaluate()

# Get optimal control
x_opt = simulator.output(0)
y_opt = simulator.output(1)
u_opt = simulator.output(2)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid,x_opt,'--')
plt.plot(tgrid,y_opt,'-')
plt.plot(tgrid,u_opt,'-.')
plt.title("Van der Pol optimization - indirect multiple shooting")
plt.xlabel('time')
plt.legend(['x trajectory','y trajectory','u trajectory'])
plt.grid()
plt.show()







