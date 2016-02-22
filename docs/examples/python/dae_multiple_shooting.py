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

"""
Solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE)

minimize     integral_{t=0}^{10} x0^2 + x1^2 + u^2  dt
x0,x1,z,u

subject to   dot(x0) == z*x0-x1+u     \
             dot(x1) == x0             }  for 0 <= t <= 10
                   0 == x1^2 + z - 1  /
             x0(t=0) == 0
             x1(t=0) == 1
             x0(t=10) == 0
             x1(t=10) == 0
             -0.75 <= u <= 1  for 0 <= t <= 10

The method used is direct multiple shooting.

Joel Andersson, 2015
"""

# Declare variables
x0 = SX.sym('x0')
x1 = SX.sym('x1')
x = vertcat(x0, x1) # Differential states
z = SX.sym('z')       # Algebraic variable
u = SX.sym('u')       # Control

# Differential equation
f_x = vertcat(z*x0-x1+u, x0)

# Algebraic equation
f_z = x1**2 + z - 1

# Lagrange cost term (quadrature)
f_q = x0**2 + x1**2 + u**2

# Create an integrator
dae = {'x':x, 'z':z, 'p':u, 'ode':f_x, 'alg':f_z, 'quad':f_q}
opts = {'tf':0.5} # interval length
I = integrator('I', 'idas', dae, opts)

# Number of intervals
nk = 20

# Start with an empty NLP
w = []   # List of variables
lbw = [] # Lower bounds on w
ubw = [] # Upper bounds on w
G = []   # Constraints
J = 0    # Cost function

# Initial conditions
Xk = MX.sym('X0', 2)
w.append(Xk)
lbw += [ 0, 1 ]
ubw += [ 0, 1 ]

# Loop over all intervals
for k in range(nk):
  # Local control
  Uk = MX.sym('U'+str(k))
  w.append(Uk)
  lbw += [-0.75]
  ubw += [ 1.00]

  # Call integrator function
  Ik = I(x0=Xk, p=Uk)
  Xk = Ik['xf']
  J = J + Ik['qf'] # Sum quadratures

  # "Lift" the variable
  X_prev = Xk
  Xk = MX.sym('X'+str(k+1), 2)
  w.append(Xk)
  lbw += [-inf, -inf]
  ubw += [ inf,  inf]
  G.append(X_prev - Xk)

# Allocate an NLP solver
nlp = {'x':vertcat(*w), 'f':J, 'g':vertcat(*G)}
opts = {'ipopt.linear_solver':'ma27'}
solver = nlpsol('solver', 'ipopt', nlp, opts)

# Pass bounds, initial guess and solve NLP
sol = solver(lbx = lbw, # Lower variable bound
             ubx = ubw,  # Upper variable bound
             lbg = 0.0,  # Lower constraint bound
             ubg = 0.0,  # Upper constraint bound
             x0  = 0.0) # Initial guess

# Plot the results
import matplotlib.pyplot as plt
plt.figure(1)
plt.clf()
plt.plot(linspace(0., 10., nk+1), sol['x'][0::3],'--')
plt.plot(linspace(0., 10., nk+1), sol['x'][1::3],'-')
plt.plot(linspace(0., 10., nk), sol['x'][2::3],'-.')
plt.title('Van der Pol optimization - multiple shooting')
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()
