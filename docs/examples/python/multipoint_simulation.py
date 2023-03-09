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
import numpy as np
import matplotlib.pyplot as plt

# Consider the following nonlinear system with two states and one control input:
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x = vertcat(x1, x2)
u = MX.sym('u')
xdot = vertcat((1-x2**2)*x1 - x2 + u, x1)

# In addition to the dynamics, we are also interested in the sum-of-squares distance from the origin:
y = x1**2 + x2**2

# We want to simulate this system over T seconds, with the control discretized into N intervals:
# Time horizon, discretization
N = 10
T = 10.
tgrid = np.linspace(0, T, N+1)

# We will integrate the problem using CVODES, while also calculating the integral of y
dae = {'x':x, 'u':u, 'ode':xdot, 'quad':y}
F = integrator('F', 'cvodes', dae, tgrid[0], tgrid[1:])

# We are also interested in the pointwise values of y, so let's create a function for that too
yfun = Function('yfun', [x], [y], ['x'], ['y'])

# Let's define initial conditions for x and values for the piecewise constant controls
u = np.linspace(-1, 1, N)
x0 = np.array([0, 0])

# We can simulate the system as follows:
Fk = F(x0 = x0, u = u)
xf = Fk['xf']
qf = Fk['qf']

# Let's add the state at the initial time to xf and qf
xf = horzcat(x0, xf)
qf = horzcat(0, qf)

# Call yfun to get y at each grid point
yf = yfun(xf)

# Convert CasADi matrices to numpy arrays
x_sim = xf.full()
y_sim = yf.full()
q_sim = qf.full()

# Plot simulation trajectories and pointwise sum of squares contribution to y, which we compare with the integral of y
fig, [ax1, ax2] = plt.subplots(nrows = 1, ncols = 2)

ax1.plot(tgrid, x_sim[0,:], '*-', label = '$x_1(t)$')
ax1.plot(tgrid, x_sim[1,:], '*-', label = '$x_2(t)$')
ax1.plot(tgrid, y_sim[0,:], '*-', label = '$(x_1(t))^2 + (x_2(t))^2$')
ax1.step(tgrid, np.insert(u, 0, np.nan), label = '$u(t)$')
ax1.set_xlabel(r'$t$')
ax1.set_title('Trajectories')
ax1.legend()
ax1.grid()

ax2.plot(tgrid, q_sim[0,:], '*-', label = r'$\int_{t=0}^{t_k}{((x_1(t))^2 + (x_2(t))^2) \, dt}$')
ax2.plot(tgrid, T/N * np.cumsum(y_sim), '*-', label = r'$\sum_{i=0}^{k}{((x_1(t_i))^2 + (x_2(t_i))^2) \, \frac{T}{N}}$')
ax2.set_xlabel(r'$t_k$')
ax2.set_title('Integral cost')
ax2.legend()
ax2.grid()



plt.show()
