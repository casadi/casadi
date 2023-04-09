# %% Build the problem (CasADi code, independent of alpaqa)
import casadi as cs
import numpy as np

# Make symbolic decision variables
x1, x2 = cs.SX.sym("x1"), cs.SX.sym("x2")
x = cs.vertcat(x1, x2)  # Collect decision variables into one vector
# Make a parameter symbol
p = cs.SX.sym("p")

# Objective function f and the constraints function g
f = (1 - x1) ** 2 + p * (x2 - x1**2) ** 2
g = cs.vertcat(
    (x1 - 0.5) ** 3 - x2 + 1,
    x1 + x2 - 1.5,
)

# Set the bounds
C = [-0.25, -0.5], [1.5, 2.5]  # -0.25 <= x1 <= 1.5, -0.5 <= x2 <= 2.5
D = [-np.inf, -np.inf], [0, 0]  #         g1 <= 0,           g2 <= 0

# %% Generate and compile C-code for the objective and constraints using alpaqa
from alpaqa import minimize

problem = (
    minimize(f, x)  #       Objective function f(x)
    .subject_to_box(C)  #   Box constraints x ∊ C
    .subject_to(g, D)  #    General ALM constraints g(x) ∊ D
    .with_param(p, [10])  # Parameter with default value (can be changed later)
).compile()

# %% Build a solver with the default parameters
import alpaqa as pa

inner_solver = pa.PANOCSolver()
solver = pa.ALMSolver(inner_solver)

# %% Build a solver with custom parameters
inner_solver = pa.PANOCSolver(
    dict(
        max_iter=50,
        stop_crit=pa.PANOCStopCrit.ApproxKKT,
        print_interval=1,
    ),
    dict(memory=10),
)

# You can attach a callback that is called on each iteration, and keeps track of
# the iterates so they can be plotted later
iterates = []
def callback(it: pa.PANOCProgressInfo):
    if it.k == 0:
        iterates.append([])
    iterates[-1].append(np.copy(it.x))
# Note: the iterate values like it.x are only valid within the callback, if you
# want to save them for later, you have to make a copy.
inner_solver.set_progress_callback(callback)

alm_params = pa.ALMParams(
    tolerance=1e-8,
    dual_tolerance=1e-8,
    initial_penalty=50,
    penalty_update_factor=20,
    print_interval=1,
)
solver = pa.ALMSolver(
    alm_params,
    inner_solver,
)

# %% Compute a solution

x0 = [0.1, 1.8]  # Initial guess for decision variables
y0 = [0.0, 0.0]  # Initial guess for Lagrange multipliers of g(x)

# Solve the problem
x_sol, y_sol, stats = solver(problem, x0, y0)

# Print the results
print(solver)
print(stats["status"])
print(f"Solution:      {x_sol}")
print(f"Multipliers:   {y_sol}")
print(f"Cost:          {problem.eval_f(x_sol)}")
from pprint import pprint
pprint(stats)

# %% Plot the results

import matplotlib.pyplot as plt
from matplotlib import patheffects

cost_function_v = np.vectorize(problem.eval_f, signature="(n)->()")
constr_function_v = np.vectorize(problem.eval_g, signature="(n)->(m)")

x = np.linspace(-1.5, 1.5, 256)
y = np.linspace(-0.5, 2.5, 256)
X, Y = np.meshgrid(x, y)
XY = np.vstack([[X], [Y]]).T

plt.figure(figsize=(10, 6))
# Draw objective function contours
Zf = cost_function_v(XY).T
plt.contourf(X, Y, Zf, 32)
plt.colorbar()
# Draw constraints
Zgc, Zgl = constr_function_v(XY).T
line_opts = dict(colors="black", linewidths=0.8, linestyles="-")
fx = [patheffects.withTickedStroke(spacing=7, linewidth=0.8)]
cgc = plt.contour(X, Y, Zgc, [0], **line_opts)
plt.setp(cgc.collections, path_effects=fx)
cgl = plt.contour(X, Y, Zgl, [0], **line_opts)
plt.setp(cgl.collections, path_effects=fx)
xl = plt.contour(X, Y, -X, [-problem.C.lowerbound[0]], **line_opts)
plt.setp(xl.collections, path_effects=fx)

plt.title("PANOC+ALM Rosenbrock example")
plt.xlabel("$x_1$")
plt.ylabel("$x_2$")

# Draw iterates and solution
for xy in map(np.array, iterates):
    plt.plot(xy[:, 0], xy[:, 1], "r:.", markersize=4, linewidth=1)
plt.plot(x_sol[0], x_sol[1], "ro", markersize=10, fillstyle="none")

plt.tight_layout()
plt.show()
