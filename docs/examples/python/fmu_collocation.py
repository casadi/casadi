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

if False:
    # Use FMU compiled from Modelica using Dymola
    dae = DaeBuilder('crane', r'../../../test/data/Crane.fmu')
else:
    # Start with an empty DaeBuilder instance
    dae = DaeBuilder('crane')

    # Hard coded constants
    mL = 1.  # Mass of the load
    mC = 1. # Mass of the crane
    g = 9.81 # Gravity
    length = 1.0 # Rope/pole length
    inverted = True # Inverted pendulum on a cart rather than crane with hanging load

    # States
    xL = dae.add('xL', dict(start = 0, description = 'Horizontal position of the load'))
    y = dae.add('y', dict(start = length if inverted else -length,
                        description = 'Vertical position of the load'))
    xC = dae.add('xC', dict(start = 0, description = 'Horizontal position of the crane'))
    dxL = dae.add('dxL', dict(start = 0, description = 'Horizontal velocity of the load'))
    dy = dae.add('dy', dict(start = 0, description = 'Vertical velocity of the load'))
    dxC = dae.add('dxC', dict(start = 0, description = 'Horizontal velocity of the crane'))

    # Input
    u = dae.add('u', 'input', dict(min = -2, max = 2, start = 0,
                                description = 'Horizontal force applied to the crane'))

    # Algebraic variables
    T = dae.add('T', dict(start = g if inverted else -g,
                        description = 'Tension in the rope/pole'))

    # Ordinary differential equations
    ddxL = (xL-xC)/length*T/mL
    ddy = g - y/length*T/mL
    ddxC = ((xL-xC)/length*T - u)/mC
    dae.eq(dae.der(xL), dxL)
    dae.eq(dae.der(y), dy)
    dae.eq(dae.der(xC), dxC)
    dae.eq(dae.der(dxL), ddxL)
    dae.eq(dae.der(dy), ddy)
    dae.eq(dae.der(dxC), ddxC)

    # Original algebraic equation
    res = y**2 + (xL-xC)**2 - length**2
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

# Simulation duration
dae.set_start_time(0)
dae.set_stop_time(5)

# DAE right-hand-side function
ffcn = dae.create('ffcn', ['x', 'z', 'u'], ['ode', 'alg', 'y'])

# Variable names
U = dae.u()
X = dae.x()
Z = dae.z()
Y = dae.y()

# Refence location for crane and load
xref = 0.1

# Quadratic cost at the beginning of control intervals
X_ref = np.array([xref if n in ('xL', 'xC') else 0 for n in X])
Y_ref = np.array([xref if n in ('xL', 'xC') else 0 for n in Y])
X_weights = np.array([1. if n in ('xL', 'xC') else 0 for n in X])
Y_weights = np.array([1. if n in ('xL', 'xC') else 0 for n in Y])

# Quadratic cost at final time
Xf_ref = np.array([xref if n in ('xL', 'xC') else 0 for n in X])
Yf_ref = np.array([xref if n in ('xL', 'xC') else 0 for n in Y])
Xf_weights = np.array([1. if n in ('xL', 'xC', 'dxL', 'dy') else 0 for n in X])
Yf_weights = np.array([1. if n in ('xL', 'xC', 'dxL', 'dy') else 0 for n in Y])

# Discretization
nk = 50 # Number of control intervals
deg = 4 # Degree of interpolating polynomial

# Size of the finite elements
h = (dae.stop_time() - dae.start_time())/nk

# We will use collocation discretization using Legendre roots, cf. 
# Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
# by Lorenz Biegler (2010).
# The roots can be queried from CasADi or looked up in the above textbook
tau_root = collocation_points(deg, 'radau')

# We can query the coefficient for the interpolating polynomials using a helper function
# C: Coefficients for calculating state derivatives
# D: Coefficients for calculating state at the end of the interval
# B: Coefficients for calculating a quadrature at the end of the interval
C, D, B = collocation_coeff(tau_root)

# We will construct an NLP incrementally, starting with no decision variables, 
# and constraints, and zero objective
W = []
G = []
J = 0

# Solution trajectories
U_all = []
X_all = []
Z_all = []
Y_all = []

# Time grids, for solution
Ugrid = []
Xgrid = []
Zgrid = []
Ygrid = []

# Current time
tk = 0

# Add a variable for the initial state
Xk = MX.sym('X0', len(X))
X_all.append(Xk)
W.append(Xk)
Xgrid.append(tk)
# Loop over all control intervals, constructing the NLP and x, u trajectories
for k in range(nk):
    # Add contribution to the cost
    Dx = Xk - X_ref
    J += dot(Dx, X_weights * Dx)
    if k > 0 and len(Y) > 0:
        # Add contribution to the cost
        Dy = Yk - Y_ref
        J += dot(Dy, Y_weights * Dy)
    # Symbolic expression for the control for the interval
    Uk = MX.sym('u' + str(k), len(U))
    U_all.append(Uk)
    W.append(Uk)
    Ugrid.append(tk)
    # Declare variables at each collocation point
    Xc = []
    Zc = []
    Yc = []
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
        # Output variables at collocation points
        Ykj = MX.sym('Y' + str(k) + '_' + str(j), len(Y))
        Yc.append(Ykj)
        Y_all.append(Ykj)
        W.append(Ykj)
        Ygrid.append(tk + tau_root[j-1]*h)
    # Collect the collocation equations
    for j in range(deg):
      # Expression for the state derivative at the collocation point
      Xdot_j = C[0,j] * Xk
      for r in range(deg): Xdot_j += C[r+1, j] * Xc[r]
      # Call the DAE right-hand-side function
      [Ode_j, Alg_j, Y_j] = ffcn(Xc[j], Zc[j], Uk)
      # Append collocation equations
      G.append(h*Ode_j - Xdot_j)
      # Add algebraic equations
      G.append(Alg_j)
      # Add output equations
      G.append(Y_j - Yc[j])
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
    # Output, algebraic variable at the end of the interval (only correct for Radau)
    Zk = Zc[-1]
    Yk = Yc[-1]

# Include end of control interval
Ugrid.append(tk)

# Add cost contribution from the final state
Dx = Xk - Xf_ref
J += dot(Dx, Xf_weights * Dx)
if len(Y) > 0:
    # Add contribution to the cost
    Dy = Yk - Yf_ref
    J += dot(Dy, Yf_weights * Dy)

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
Y_all = hcat(Y_all)
to_W = Function('to_W', [X_all, Z_all, Y_all, U_all], [W], ['X', 'Z', 'Y', 'U'], ['W'])
from_W = Function('from_W', [W], [X_all, Z_all, Y_all, U_all], ['W'], ['X', 'Z', 'Y', 'U'])

# Variable bounds, initial guess: Get data structures of correct size
[X0, Z0, Y0, U0] = from_W(0)
[lbX, lbZ, lbY, lbU] = from_W(-np.inf)
[ubX, ubZ, ubY, ubU] = from_W(np.inf)

# Initial guess, bounds for control trajectory
Uinit = dae.get(U)
Umin = dae.min(U)
Umax = dae.max(U)
for i, v0 in enumerate(Uinit):
    U0[i, :] = v0
    lbU[i, :] = Umin[i]
    ubU[i, :] = Umax[i]

# Initial guess, bounds for state trajectory
Xinit = dae.get(X)
Xmin = dae.min(X)
Xmax = dae.max(X)
for i, v0 in enumerate(Xinit):
    X0[i, :] = v0
    lbX[i, :] = Xmin[i]
    lbX[i, 0] = v0
    ubX[i, :] = Xmax[i]
    ubX[i, 0] = v0

# Initial guess for algebraic variable trajectory
Zinit = dae.get(Z)
Zmin = dae.min(Z)
Zmax = dae.max(Z)
for i, v0 in enumerate(Zinit):
    Z0[i, :] = v0
    lbZ[i, :] = Zmin[i]
    ubZ[i, :] = Zmax[i]

# Initial guess for output trajectory
Yinit = dae.get(Y)
Ymin = dae.min(Y)
Ymax = dae.max(Y)
for i, v0 in enumerate(Yinit):
    Y0[i, :] = v0
    lbY[i, :] = Ymin[i]
    ubY[i, :] = Ymax[i]

# Solve the NLP
W0 = to_W(X0, Z0, Y0, U0)
lbW = to_W(lbX, lbZ, lbY, lbU)
ubW = to_W(ubX, ubZ, ubY, ubU)
sol = solver(x0 = W0, lbx = lbW, ubx = ubW, lbg = 0, ubg = 0)
Xopt, Zopt, Yopt, Uopt = from_W(sol['x'])
Xopt = Xopt.full()
Zopt = Zopt.full()
Yopt = Yopt.full()
Uopt = Uopt.full()

# Variables to plot
plot_vars = X + Z + Y + U

# Create figure and axes
n_fig = len(plot_vars)
ncols = min(4, int(np.ceil(np.sqrt(n_fig))))  # up to 4 side-by-side
nrows = int(np.ceil(n_fig / ncols))
figwidth = 20
figheight = max(10, nrows*4)
fig = plt.figure(figsize = (figwidth, figheight), layout = 'constrained')
ax = []
for i in range(n_fig):
    ax.append(fig.add_subplot(nrows, ncols, i + 1))

# Plot trajectories
for k, n in enumerate(plot_vars):
    # Allow initial underscore in names, if any
    name = r"$\_$" + n[1:] if n.startswith('_') else n
    # Plot trajectories
    if n in X:
        i = X.index(n)
        ax[k].plot(Xgrid, Xopt[i, :], label=f'x[{i}]: ' + name)
    elif n in Z:
        i = Z.index(n)
        ax[k].plot(Zgrid, Zopt[i, :], label=f'z[{i}]: ' + name)
    elif n in Y:
        i = Y.index(n)
        ax[k].plot(Ygrid, Yopt[i, :], label=f'y[{i}]: ' + name)
    else:
        assert n in U
        i = U.index(n)
        ax[k].step(Ugrid, np.insert(Uopt[i, :], 0, np.nan), label=f'u[{i}]: ' + name)
    ax[k].legend()
    ax[k].grid()
plt.show()
