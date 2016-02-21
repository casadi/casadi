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

nk = 20    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
u  = SX.sym('u')    # control
x  = SX.sym('x',2)  # states

# ODE right hand side and quadratures
xdot = vertcat((1 - x[1]*x[1])*x[0] - x[1] + u, x[0])
qdot = x[0]*x[0] + x[1]*x[1] + u*u

# Create an integrator
dae = {'x':x, 'p':u, 'ode':xdot, 'quad':qdot}
F = integrator('F', 'cvodes', dae, {'tf':tf/nk})

# All controls (use matrix graph)
U = [MX.sym('u' + str(k)) for k in range(nk)]

# The initial state (x_0=0, x_1=1)
x0 = [0,1]
X  = MX(x0)

# Objective function
f = 0

# Build a graph of integrator calls
for k in range(nk):
  res = F(x0=X, p=U[k])
  X = res['xf']
  f += res['qf']

# Terminal constraints: x_0(T)=x_1(T)=0
g = X

# Allocate an NLP solver
nlp = {'x':vertcat(*U), 'f':f, 'g':g}
solver = nlpsol('solver', 'ipopt', nlp)

# Solve the problem
sol = solver(lbx = -0.75,
             ubx = 1,
             x0 = 0,
             lbg = 0,
             ubg = 0)
              
# Retrieve the solution
u_opt = sol['x']

# Simulate to get optimal X trajectory
x_opt = [x0]
for k in range(nk):
  res = F(x0=x_opt[-1], p=u_opt[k])
  x_opt.append(res['xf'])
x0_opt = [xk[0] for xk in x_opt]
x1_opt = [xk[1] for xk in x_opt]

# Time grid for printing
tgrid = [tf/nk*k for k in range(nk+1)]

# Plot the results
import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x0_opt, '--')
plt.plot(tgrid, x1_opt, '-')
plt.step(tgrid, vertcat(DM.nan(1), u_opt), '-.') # Note: first entry is ignored
plt.title('Van der Pol optimization - single shooting')
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

