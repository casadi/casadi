## @example simple_optimization/python/rosenbrock.py
# This code contains a minimal example of an optimization problem that can be 
# built and solved using `alpaqa`, it includes visualization of the iterates.

# %% Build the problem for PANOC+ALM (CasADi code, independent of alpaqa)
import casadi as cs

# Make symbolic decision variables
x1, x2 = cs.SX.sym("x1"), cs.SX.sym("x2")
# Make a parameter symbol
p = cs.SX.sym("p")

# Expressions for the objective function f and the constraints g
f_expr = (1 - x1) ** 2 + p * (x2 - x1 ** 2) ** 2
g_expr = cs.vertcat(
    (x1 - 0.5) ** 3 - x2 + 1,
    x1 + x2 - 1.5,
)

# Collect decision variables into one vector
x = cs.vertcat(x1, x2)
# Convert the symbolic expressions to CasADi functions
f = cs.Function("f", [x, p], [f_expr])
g = cs.Function("g", [x, p], [g_expr])

# %% Generate and compile C-code for the objective and constraints using alpaqa
import alpaqa as pa

# Compile and load the problem
prob = pa.generate_and_compile_casadi_problem(f, g)

# Set the bounds
import numpy as np
prob.C.lowerbound = [-0.25, -0.5]       # -0.25 <= x1 <= 1.5
prob.C.upperbound = [1.5, 2.5]          # -0.5  <= x2 <= 2.5
prob.D.lowerbound = [-np.inf, -np.inf]  # g1 <= 0
prob.D.upperbound = [0, 0]              # g2 <= 0

# Set parameter to some value
prob.param = [10.]

# %% Build an inner Structured PANOC solver with custom parameters
innersolver = pa.StructuredPANOCLBFGSSolver(
    pa.StructuredPANOCLBFGSParams(
        max_iter=1000,
        stop_crit=pa.PANOCStopCrit.ApproxKKT,
    ),
    pa.LBFGSParams(
        memory=10,
    ),
)

# You can attach a callback that is called on each iteration, and keeps track of
# the iterates so they can be plotted later
iterates = []
def cb(it): iterates.append(np.copy(it.x))
# Note: the iterate values like it.x are only valid within the callback, if you
# want to save them for later, you have to make a copy.
innersolver.set_progress_callback(cb)

# %% Make an ALM solver with default parameters, using the PANOC solver
solver = pa.ALMSolver(
    pa.ALMParams(
        ε=1e-10,
        δ=1e-10,
        Σ_0=0,
        σ_0=2,
        Δ=20,
    ),
    innersolver
)

# %% Compute a solution 

# Set initial guesses at arbitrary values
x0 = np.array([0.1, 1.8]) # decision variables
y0 = np.zeros((prob.m,))  # Lagrange multipliers for g(x)

# Solve the problem
x_sol, y_sol, stats = solver(prob, x0, y0)

# Print the results
print(stats["status"])
print(f"Solution:      {x_sol}")
print(f"Multipliers:   {y_sol}")
print(f"Cost:          {prob.f(x_sol)}")
from pprint import pprint
pprint(stats)

# %% Plot the results

import matplotlib.pyplot as plt
from matplotlib import patheffects

cost_function_v = np.vectorize(prob.f, signature='(n)->()')
constraint_g_v = np.vectorize(prob.g, signature='(n)->(m)')

x = np.linspace(-1.5, 1.5, 256)
y = np.linspace(-0.5, 2.5, 256)
X, Y = np.meshgrid(x, y)
XY = np.vstack([[X], [Y]]).T

plt.figure(figsize=(10, 6))
# Draw objective function
Zf = cost_function_v(XY).T
plt.contourf(X, Y, Zf, 32)
plt.colorbar()
# Draw constraints
Zg = constraint_g_v(XY)
Zgc = Zg[:,:,0].T
Zgl = Zg[:,:,1].T
fx = [patheffects.withTickedStroke(spacing=7, linewidth=0.8)]
cgc = plt.contour(X, Y, Zgc, [0], colors='black', linewidths=0.8, linestyles='-')
plt.setp(cgc.collections, path_effects=fx)
cgl = plt.contour(X, Y, Zgl, [0], colors='black', linewidths=0.8, linestyles='-')
plt.setp(cgl.collections, path_effects=fx)
xl = plt.contour(X, Y, -X, [-prob.C.lowerbound[0]], colors='black', linewidths=0.8, linestyles='-')
plt.setp(xl.collections, path_effects=fx)

plt.title("PANOC+ALM Rosenbrock example")
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")

# Draw iterates and solution
xy = np.array(iterates)
plt.plot(xy[:,0], xy[:,1], 'r:.', markersize=4, linewidth=1)
plt.plot(x_sol[0], x_sol[1], 'ro', markersize=10, fillstyle='none')

plt.tight_layout()
plt.show()
