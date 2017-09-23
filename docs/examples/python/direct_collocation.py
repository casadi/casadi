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
import casadi as c
import numpy as n
import matplotlib.pyplot as plt

# Degree of interpolating polynomial
d = 3

# Get collocation points
tau_root = n.append(0, c.collocation_points(d, 'legendre'))

# Coefficients of the collocation equation
C = n.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = n.zeros(d+1)

# Coefficients of the quadrature function
B = n.zeros(d+1)

# Construct polynomial basis
for j in range(d+1):
    # Construct Lagrange polynomials to get the polynomial basis at the collocation point
    p = n.poly1d([1])
    for r in range(d+1):
        if r != j:
            p *= n.poly1d([1, -tau_root[r]]) / (tau_root[j]-tau_root[r])

    # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D[j] = p(1.0)

    # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = n.polyder(p)
    for r in range(d+1):
        C[j,r] = pder(tau_root[r])

    # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = n.polyint(p)
    B[j] = pint(1.0)

# Time horizon
T = 10.

# Declare model variables
x1 = c.SX.sym('x1')
x2 = c.SX.sym('x2')
x = c.vertcat(x1, x2)
u = c.SX.sym('u')

# Model equations
xdot = c.vertcat((1-x2**2)*x1 - x2 + u, x1)

# Objective term
L = x1**2 + x2**2 + u**2

# Continuous time dynamics
f = c.Function('f', [x, u], [xdot, L], ['x', 'u'], ['xdot', 'L'])

# Control discretization
N = 20 # number of control intervals
h = T/N

# Start with an empty NLP
w=[]
w0 = []
lbw = []
ubw = []
J = 0
g=[]
lbg = []
ubg = []

# For plotting x and u given w
x_plot = []
u_plot = []

# "Lift" initial conditions
X0 = c.MX.sym('X0', 2)
w.append(X0)
lbw.extend([0, 1])
ubw.extend([0, 1])
w0.extend([0, 1])
x_plot.append(X0)

# Formulate the NLP
Xk = c.MX([0, 1])
for k in range(N):
    # New NLP variable for the control
    Uk = c.MX.sym('U_' + str(k))
    w.append(Uk)
    lbw.extend([-1])
    ubw.extend([1])
    w0.extend([0])
    u_plot.append(Uk)

    # State at collocation points
    Xc = []
    for j in range(d):
        Xkj = c.MX.sym('X_'+str(k)+'_'+str(j), 2)
        Xc.append(Xkj)
        w.append(Xkj)
        lbw.extend([-0.25, -n.inf])
        ubw.extend([n.inf,  n.inf])
        w0.extend([0, 0])

    # Loop over collocation points
    Xk_end = D[0]*Xk
    for j in range(1,d+1):
       # Expression for the state derivative at the collocation point
       xp = C[0,j]*Xk
       for r in range(d): xp = xp + C[r+1,j]*Xc[r]

       # Append collocation equations
       fj, qj = f(Xc[j-1],Uk)
       g.append(h*fj - xp)
       lbg.extend([0, 0])
       ubg.extend([0, 0])

       # Add contribution to the end state
       Xk_end = Xk_end + D[j]*Xc[j-1];

       # Add contribution to quadrature function
       J = J + B[j]*qj*h

    # New NLP variable for state at end of interval
    Xk = c.MX.sym('X_' + str(k+1), 2)
    w.append(Xk)
    lbw.extend([-0.25, -n.inf])
    ubw.extend([n.inf,  n.inf])
    w0.extend([0, 0])
    x_plot.append(Xk)

    # Add equality constraint
    g.append(Xk_end-Xk)
    lbg.extend([0, 0])
    ubg.extend([0, 0])

# Create an NLP solver
prob = {'f': J, 'x': c.vertcat(*w), 'g': c.vertcat(*g)}
solver = c.nlpsol('solver', 'ipopt', prob);

# Function to get x and u trajectories from w
trajectories = c.Function('trajectories',
     [c.vertcat(*w)], [c.horzcat(*x_plot), c.horzcat(*u_plot)], ['w'], ['x', 'u'])

# Solve the NLP
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
x_opt, u_opt = trajectories(sol['x'])
x_opt = x_opt.full() # to numpy array
u_opt = u_opt.full() # to numpy array

# Plot the result
tgrid = n.linspace(0, T, N+1)
plt.figure(1)
plt.clf()
plt.plot(tgrid, x_opt[0], '--')
plt.plot(tgrid, x_opt[1], '-')
plt.step(tgrid, n.append(n.nan, u_opt[0]), '-.')
plt.xlabel('t')
plt.legend(['x1','x2','u'])
plt.grid()
plt.show()
