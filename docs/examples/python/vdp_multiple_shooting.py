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
coll = False # Use collocation integrator

# Declare variables (use scalar graph)
t  = SX.sym("t")    # time
u  = SX.sym("u")    # control
x  = SX.sym("x",3)  # state

# ODE rhs function
ode = vertcat([(1 - x[1]*x[1])*x[0] - x[1] + u, \
       x[0], \
       x[0]*x[0] + x[1]*x[1] + u*u])
dae = SXFunction(daeIn(x=x,p=u,t=t),daeOut(ode=ode))

# Create an integrator
if coll:
  integrator = OldCollocationIntegrator(dae)
  integrator.setOption("number_of_finite_elements",5)
  integrator.setOption("interpolation_order",5)
  integrator.setOption("collocation_scheme","legendre")
  integrator.setOption("implicit_solver","kinsol")
  integrator.setOption("implicit_solver_options",\
    {'linear_solver' : CSparse})
  integrator.setOption("expand_f",True)
else:
  integrator = Integrator("cvodes", dae)
  integrator.setOption("abstol",1e-8) # tolerance
  integrator.setOption("reltol",1e-8) # tolerance
  integrator.setOption("steps_per_checkpoint",1000)

integrator.setOption("tf",tf/nk) # final time
integrator.init()

# Total number of variables
nv = 1*nk + 3*(nk+1)

# Declare variable vector
V = MX.sym("V", nv)

# Get the expressions for local variables
U = V[0:nk]
X0 = V[nk+0:nv:3]
X1 = V[nk+1:nv:3]
X2 = V[nk+2:nv:3]

# Variable bounds initialized to +/- inf
VMIN = -inf*NP.ones(nv)
VMAX = inf*NP.ones(nv)

# Control bounds
VMIN[0:nk] = -0.75
VMAX[0:nk] = 1.0

# Initial condition
VMIN[nk+0] = VMAX[nk+0] = 0
VMIN[nk+1] = VMAX[nk+1] = 1
VMIN[nk+2] = VMAX[nk+2] = 0

# Terminal constraint
VMIN[nv-3] = VMAX[nv-3] = 0
VMIN[nv-2] = VMAX[nv-2] = 0

# Initial solution guess
VINIT = NP.zeros(nv)

# Constraint function with bounds
g = [];  g_min = []; g_max = []

# Build up a graph of integrator calls
for k in range(nk):
  # Local state vector
  Xk = vertcat((X0[k],X1[k],X2[k]))
  Xk_next = vertcat((X0[k+1],X1[k+1],X2[k+1]))
  
  # Call the integrator
  Xk_end, = integratorOut(integrator(integratorIn(x0=Xk,p=U[k])),"xf")
  
  # append continuity constraints
  g.append(Xk_next - Xk_end)
  g_min.append(NP.zeros(Xk.size()))
  g_max.append(NP.zeros(Xk.size()))

# Objective function: L(T)
f = X2[nk]

# Continuity constraints: 0<= x(T(k+1)) - X(T(k)) <=0
g = vertcat(g)

# Create NLP solver instance
nlp = MXFunction(nlpIn(x=V),nlpOut(f=f,g=g))
solver = NlpSolver("ipopt", nlp)

#solver.setOption("verbose",True)
solver.init()

# Set bounds and initial guess
solver.setInput(VMIN,  "lbx")
solver.setInput(VMAX,  "ubx")
solver.setInput(VINIT, "x0")
solver.setInput(NP.concatenate(g_min),"lbg")
solver.setInput(NP.concatenate(g_max),"ubg")

# Solve the problem
solver.evaluate()

# Retrieve the solution
v_opt = solver.getOutput("x")
u_opt = v_opt[0:nk]
x0_opt = v_opt[nk+0::3]
x1_opt = v_opt[nk+1::3]
x2_opt = v_opt[nk+2::3]

# Get values at the beginning of each finite element
tgrid_x = NP.linspace(0,10,nk+1)
tgrid_u = NP.linspace(0,10,nk)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid_x,x0_opt,'--')
plt.plot(tgrid_x,x1_opt,'-')
plt.plot(tgrid_u,u_opt,'-.')
plt.title("Van der Pol optimization - multiple shooting")
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

