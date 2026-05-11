#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# -*- coding: utf-8 -*-
"""
Example on how to formulate and solve an optimal control problem using direct collocation.
The model is a semi-explicit DAE.

Based on collocation example by Mario Zanon, Sebastien Gross 2012
Reformulated to a semi-explict formulation and modern syntax by Joel Andersson, 2026
"""
from casadi import *
import numpy as np
import matplotlib.pyplot as plt

# Dynamic model: 
m = 1.
M = 1.
g = 9.81

# Control
u  = SX.sym("u")
Uv = u
U = ['u']

# Differential states
x = SX.sym("x")
y = SX.sym("y")
w = SX.sym("w")
dx = SX.sym("dx")
dy = SX.sym("dy")
dw = SX.sym("dw")
Xv = vertcat(x,y,w,dx,dy,dw)
X = ['x', 'y', 'w', 'dx', 'dy', 'dw']

# Algebraic variables
lam = SX.sym("lam")
Zv = lam
Z = ['lam']

# Ordinary differential equations
ddx = (x-w)*lam/m
ddy = g - y*lam/m
ddw = ((x-w)*lam - u)/M
Xdot = vertcat(dx, dy, dw, ddx, ddy, ddw)

# Algebraic equation
Alg = (x-w)*(ddx-ddw) + y*ddy + dy*dy + (dx-dw)*(dx-dw)

# Quadrature
xref = 0.1 # chariot reference
lagrange_term = (x-xref)**2 + (w-xref)**2
Quad = lagrange_term

# Output
mayer_term = (x-xref)*(x-xref) + (w-xref)*(w-xref) + dx*dx + dy*dy
Out = mayer_term

# DAE rhs
ffcn = Function('ffcn', [Xv, Zv, Uv], [Xdot, Alg, Quad])

# Output function
hfcn = Function('hfun', [Xv, Uv], [Out])

# -----------------------------------------------------------------------------
# Collocation setup
# -----------------------------------------------------------------------------
yl = 1. #- -> crane, + -> pendulum
tf = 5.0
nk = 50

# Degree of interpolating polynomial
deg = 4

# We will use collocation discretization using Legendre roots, cf. 
# Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
# by Lorenz Biegler (2010).
# The roots can be queried from CasADi or looked up in the above textbook
tau_root = collocation_points(deg, 'radau')

# We can query the coefficient for the interpolating polynomials using a helper function
C, D, B = collocation_coeff(tau_root)

# Size of the finite elements
h = tf/nk

# We will construct an NLP incrementally, starting with no decision variables, 
# and constraints, and zero objective
W = []
G = []
J = 0

# Solution trajectories
U_all = []
X_all = []
Z_all = []

# Time grids, for solution
Ugrid = []
Xgrid = []
Zgrid = []

# Current time
tk = 0

# Add a variable for the initial state
Xk = MX.sym('X0', len(X))
X_all.append(Xk)
W.append(Xk)
Xgrid.append(tk)
# Loop over all control intervals, constructing the NLP and x, u trajectories
for k in range(nk):
    # Symbolic expression for the control for the interval
    Uk = MX.sym('u' + str(k), len(U))
    U_all.append(Uk)
    W.append(Uk)
    Ugrid.append(tk)
    # Declare variables at each collocation point
    Xc = []
    Zc = []
    for j in range(1, deg + 1):
        # Differential state at collocation points
        Xkj = MX.sym('X' + str(k) + '_' + str(j), len(X))
        Xc.append(Xkj)
        X_all.append(Xkj)
        W.append(Xkj)
        Xgrid.append(tk + tau_root[j-1]*h)
        # Algebraic variables at collocation points
        Zkj = MX.sym('Z' + str(k) + '_' + str(j), len(Z))
        Zc.append(Zkj)
        Z_all.append(Zkj)
        W.append(Zkj)
        Zgrid.append(tk + tau_root[j-1]*h)
    # Collect the collocation equations
    for j in range(deg):
      # Expression for the state derivative at the collocation point
      Xdot_j = C[0,j] * Xk
      for r in range(deg): Xdot_j += C[r+1, j] * Xc[r]
      # Call the DAE right-hand-side function
      [Ode_j, Alg_j, Q_j] = ffcn(Xc[j], Zc[j], Uk)
      # Append collocation equations
      G.append(h*Ode_j - Xdot_j)
      # Add algebraic equations
      G.append(Alg_j)
      # Add Lagrange term contribution
      J += h * B[j] * Q_j
    # State at the end of the interval
    Xk_end = D[0] * Xk
    for j in range(deg): Xk_end = Xk_end + D[j+1] * Xc[j]
    # New variable for the state at the end of the interval
    tk = tk + h
    Xk = MX.sym('X' + str(k + 1), len(X))
    X_all.append(Xk)
    W.append(Xk)
    Xgrid.append(tk)
    # Enforce continuity
    G.append(Xk_end - Xk)

# Include end of control interval
Ugrid.append(tk)

# Add Mayer term contribution
J += hfcn(Xk, Uk)

# Concatenate vectors
W = vcat(W)
G = vcat(G)

# Form the NLP
nlp = dict(x = W, f = J, g = G)

# NLP solver options
opts = {}
opts["ipopt.max_iter"] = 100
opts["ipopt.tol"] = 1e-4
#opts["ipopt.linear_solver"] = 'ma27'

# Allocate an NLP solver
solver = nlpsol("solver", "ipopt", nlp, opts)

# Create mappings from (X, Z, U) -> (W) and back
X_all = hcat(X_all)
Z_all = hcat(Z_all)
U_all = hcat(U_all)
to_W = Function('to_W', [X_all, Z_all, U_all], [W], ['X', 'Z', 'U'], ['W'])
from_W = Function('from_W', [W], [X_all, Z_all, U_all], ['W'], ['X', 'Z', 'U'])

# Variable bounds, initial guess
[lbX, lbZ, lbU] = from_W(-np.inf)
[ubX, ubZ, ubU] = from_W(np.inf)
[X0, Z0, U0] = from_W(0)

# Specify control bounds
lbU[U.index('u'), :] = -2
ubU[U.index('u'), :] = 2

# Fixed initial condition on state
lbX[:, 0] = 0
ubX[:, 0] = 0
lbX[X.index('y'), 0] = yl
ubX[X.index('y'), 0] = yl

# Initialize y to length
X0[X.index('y'), :] = yl

# Initialize lam to gravity
Z0[Z.index('lam'), :] = sign(yl) * 9.81

# Solve the NLP
W0 = to_W(X0, Z0, U0)
lbW = to_W(lbX, lbZ, lbU)
ubW = to_W(ubX, ubZ, ubU)
sol = solver(x0 = W0, lbx = lbW, ubx = ubW, lbg = 0, ubg = 0)
Xopt, Zopt, Uopt = from_W(sol['x'])
Xopt = Xopt.full()
Zopt = Zopt.full()
Uopt = Uopt.full()

# Solution figure
fig, ax = plt.subplots(4, 2, sharex=True)
ax = ax.flatten()

# Plot trajectories
for k, n in enumerate(X + Z + U):
    # Plot trajectories
    if n in X:
        ax[k].plot(Xgrid, Xopt[X.index(n), :], label=n)
    elif n in Z:
        ax[k].plot(Zgrid, Zopt[Z.index(n), :], label=n)
    else:
        assert n in U
        ax[k].step(Ugrid, np.insert(Uopt[U.index(n), :], 0, np.nan), label=n)
    ax[k].legend()
    ax[k].grid()
plt.show()
