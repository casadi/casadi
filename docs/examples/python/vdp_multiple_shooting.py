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

nk = 20      # Control discretization
tf = 10.0    # End time

# Declare variables (use scalar graph)
u  = SX.sym('u')    # control
x  = SX.sym('x',2)  # state

# ODE rhs function and quadratures
xdot = vertcat((1 - x[1]*x[1])*x[0] - x[1] + u, x[0])
qdot = x[0]*x[0] + x[1]*x[1] + u*u

# Create an integrator
dae = {'x':x, 'p':u, 'ode':xdot, 'quad':qdot}
opts = {'tf':tf/nk} # final time
F = integrator('F', 'cvodes', dae, opts)

# Path constraints
u_min = [-0.75]
u_max = [ 1.00]

# Path constraints
x_min = [-inf, -inf]
x_max = [ inf,  inf]

# Initial and terminal constraints
x0_min = [0., 1.]
x0_max = [0., 1.]
xf_min = [0., 0.]
xf_max = [0., 0.]

# Initial guess
u_init = [0.]
x_init = [0., 0.]

# All local controls
U = [MX.sym('u' + str(k)) for k in range(nk)]

# State at shooting nodes
X = [MX.sym('x' + str(k), 2) for k in range(nk+1)]

# NLP variables with bounds and initial guess
w = []; w_min = []; w_max = []; w_init = []

# NLP constraints with bounds
g = [];  g_min = []; g_max = []

# NLP objective
J = 0.

# Treat initial state as a decision variable
w += [X[0]]
w_min += x0_min
w_max += x0_max
w_init += x_init

# Loop over control intervals
for k in range(nk):
  # Add local control to NLP variables
  w.append(U[k])
  w_min += u_min
  w_max += u_max
  w_init += u_init

  # Add state at end of control interval to NLP
  w.append(X[k+1])
  if k==nk:
    w_min += xf_min
    w_max += xf_max
  else:
    w_min += x_min
    w_max += x_max
  w_init += x_init

  # Add continuity constraints to NLP
  Fk = F(x0=X[k], p=U[k])
  g.append(X[k+1] - Fk['xf'])
  g_min += [0., 0.]
  g_max += [0., 0.]

  # Add contribution to the objective
  J += Fk['qf']

# Formulate the NLP
nlp = {'x':vertcat(*w), 'f':J, 'g':vertcat(*g)}

# Create NLP solver instance
solver = nlpsol('solver', 'ipopt', nlp)

# Solve the problem
sol = solver(lbx = w_min,
             ubx = w_max,
             x0  = w_init,
             lbg = g_min,
             ubg = g_max)

# Retrieve the solution
w_opt = sol['x']
x0_opt = w_opt[0::2+1]
x1_opt = w_opt[1::2+1]
u_opt = w_opt[2::2+1]

# Time grid for printing
tgrid = [tf/nk*k for k in range(nk+1)]

# Plot the results
import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x0_opt, '--')
plt.plot(tgrid, x1_opt, '-')
plt.step(tgrid, vertcat(DM.nan(1), u_opt), '-.') # Note: first entry is ignored
plt.title('Van der Pol optimization - multiple shooting')
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

