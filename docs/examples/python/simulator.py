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

# Time horizon, discretization
N = 50
T = 10.

# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x = vertcat(x1, x2)
u = MX.sym('u')

# Model equations
xdot = vertcat((1-x2**2)*x1 - x2 + u, x1)

# Objective term
L = x1**2 + x2**2 + u**2

# Time grid
tgrid = [T/N*k for k in range(N+1)]

# CVODES from the SUNDIALS suite
dae = {'x':x, 'u':u, 'ode':xdot, 'ydef' : L}
F = simulator('F', 'cvodes', dae, tgrid, dict(newton_scheme = 'direct', print_stats = True))

# Test input
import numpy as np
u = np.linspace(-1, 1, N)

# Simulate
Fk = F(x0 = 0, u = u)
x_sim = Fk['x'].full()
y_sim = Fk['y'].full()

import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(tgrid, x_sim[0,:], label = 'x1')
plt.plot(tgrid, x_sim[1,:], label = 'x2')
plt.plot(tgrid, y_sim[0,:], label = 'L')
plt.step(tgrid, np.insert(u, 0, np.nan), label = 'u')
plt.xlabel('t')
plt.legend()
plt.grid()
plt.show()
