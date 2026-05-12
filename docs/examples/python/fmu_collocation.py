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

# Start with an empty DaeBuilder instance
dae = DaeBuilder('crane', '', dict(detect_quad = True))

# Hard coded constants
m = 1.  # Mass of the load
M = 1. # Mass of the crane
g = 9.81
xref = 0.1 # chariot reference
length = 1.0 # Rope/pole length
inverted = True # Makes it an inverted pendulum on a cart
y0 = length if inverted else -length

# States
x = dae.add('x', dict(start = 0))
y = dae.add('y', dict(start = y0))
w = dae.add('w', dict(start = 0))
dx = dae.add('dx', dict(start = 0))
dy = dae.add('dy', dict(start = 0))
dw = dae.add('dw', dict(start = 0))

# Input
u = dae.add('u', 'input', dict(min = -2, max = 2, start = 0))

# Algebraic variables
T = dae.add('T', dict(start = g if inverted else -g))
# Ordinary differential equations
ddx = (x-w)/length*T/m
ddy = g - y/length*T/m
ddw = ((x-w)/length*T - u)/M
dae.eq(dae.der(x), dx)
dae.eq(dae.der(y), dy)
dae.eq(dae.der(w), dw)
dae.eq(dae.der(dx), ddx)
dae.eq(dae.der(dy), ddy)
dae.eq(dae.der(dw), ddw)

# Quadrature
lagrange_term = dae.add('lagrange_term')
dae.eq(dae.der(lagrange_term), (x-xref)**2 + (w-xref)**2)

# Output
mayer_term = dae.add('mayer_term', 'output')
dae.eq(mayer_term, (x-xref)*(x-xref) + (w-xref)*(w-xref) + dx*dx + dy*dy)

# Original algebraic equation
res = y**2 + (x-w)**2 - length**2
print(f'Original residual equation: {res}')

# Index reduction: Differentiate with respect to time
# In CasADi, the most efficient way to do index reduction it to use forward
# mode AD, e.g. us using Jacobian-times-vector products (jtimes):
# Using the chain rule we have: d(res)/dt = d(res)/dx * dx/dt
res1 = jtimes(res, vcat(dae.inputs('x')), vcat(dae.outputs('ode')))
res2 = jtimes(res1, vcat(dae.inputs('x')), vcat(dae.outputs('ode')))
print(f'Residual equation after index reduction: {res2}')

# Add algebraic equations
dae.eq(0, res2)

# DAE rhs
ffcn = dae.create('ffcn', ['x', 'z', 'u'], ['ode', 'alg', 'quad'])

# Output function
hfcn = dae.create('hfun', ['x', 'u'], ['y'])

# Variable names
U = dae.u()
X = dae.x()
Z = dae.z()

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
tf = 5.0
nk = 50
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

# Allocate an NLP solver
nlp = dict(x = W, f = J, g = G)
opts = {}
opts["ipopt.max_iter"] = 100
opts["ipopt.tol"] = 1e-4
solver = nlpsol("solver", "ipopt", nlp, opts)

# Create mappings from (X, Z, U) -> (W) and back
X_all = hcat(X_all)
Z_all = hcat(Z_all)
U_all = hcat(U_all)
to_W = Function('to_W', [X_all, Z_all, U_all], [W], ['X', 'Z', 'U'], ['W'])
from_W = Function('from_W', [W], [X_all, Z_all, U_all], ['W'], ['X', 'Z', 'U'])

# Variable bounds, initial guess: Get data structures of correct size
[X0, Z0, U0] = from_W(0)
[lbX, lbZ, lbU] = from_W(-np.inf)
[ubX, ubZ, ubU] = from_W(np.inf)

# Initial guess, bounds for control trajectory
Uinit = dae.start(U)
Umin = dae.min(U)
Umax = dae.max(U)
for i, v0 in enumerate(Uinit):
    U0[i, :] = v0
    lbU[i, :] = Umin[i]
    ubU[i, :] = Umax[i]

# Initial guess, bounds for state trajectory
Xinit = dae.start(X)
Xmin = dae.min(X)
Xmax = dae.max(X)
for i, v0 in enumerate(Xinit):
    X0[i, :] = v0
    lbX[i, :] = Xmin[i]
    lbX[i, 0] = v0
    ubX[i, :] = Xmax[i]
    ubX[i, 0] = v0

# Initial guess for algebraic variable trajectory
Zinit = dae.start(Z)
Zmin = dae.min(Z)
Zmax = dae.max(Z)
for i, v0 in enumerate(Zinit):
    Z0[i, :] = v0
    lbZ[i, :] = Zmin[i]
    ubZ[i, :] = Zmax[i]

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
