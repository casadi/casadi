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
x  = SX.sym('x',3)  # state

# ODE rhs function
ode = vertcat([(1 - x[1]*x[1])*x[0] - x[1] + u, \
       x[0], \
       x[0]*x[0] + x[1]*x[1] + u*u])
dae = {'x':x, 'p':u, 'ode':ode}

# Create an integrator
opts = {'tf':tf/nk} # final time
opts['abstol'] = 1e-8 # tolerance
opts['reltol'] = 1e-8 # tolerance
opts['steps_per_checkpoint'] = 1000
F = integrator('F', 'cvodes', dae, opts)

# Path constraints
u_min = [-0.75]
u_max = [ 1.00]

# Path constraints
x_min = [-inf, -inf, -inf]
x_max = [ inf,  inf,  inf]

# Initial and terminal constraints
x0_min = [0., 1., 0.]
x0_max = [0., 1., 0.]
xf_min = [0., 0., -inf]
xf_max = [0., 0.,  inf]

# Initial guess
u_init = [0.]
x_init = [0., 0., 0.]

# All local controls
U = [MX.sym('u' + str(k)) for k in range(nk)]

# State at shooting nodes
X = [MX.sym('x' + str(k), 3) for k in range(nk+1)]

# NLP variables with bounds and initial guess
w = []; w_min = []; w_max = []; w_init = []

# NLP constraints with bounds
g = [];  g_min = []; g_max = []

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
  Fk = F({'x0':X[k],'p':U[k]})
  g.append(X[k+1] - Fk['xf'])
  g_min += [0., 0., 0.]
  g_max += [0., 0., 0.]

# Formulate the NLP
nlp = {'x':vertcat(w), 'f':X[-1][2], 'g':vertcat(g)}

# Create NLP solver instance
solver = nlpsol('solver', 'ipopt', nlp)

# Solve the problem
sol = solver({'lbx' : w_min,
              'ubx' : w_max,
              'x0' : w_init,
              'lbg' : g_min,
              'ubg' : g_max})

# Retrieve the solution
w_opt = sol['x']
x0_opt = w_opt[0::3+1]
x1_opt = w_opt[1::3+1]
x2_opt = w_opt[2::3+1]
u_opt = w_opt[3::3+1]

# Time grid for printing
import numpy
tgrid_x = numpy.linspace(0, 10, nk+1)
tgrid_u = numpy.linspace(0, 10, nk)

# Plot the results
import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid_x, x0_opt, '--')
plt.plot(tgrid_x, x1_opt, '-')
plt.plot(tgrid_u, u_opt, '-.')
plt.title('Van der Pol optimization - multiple shooting')
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

