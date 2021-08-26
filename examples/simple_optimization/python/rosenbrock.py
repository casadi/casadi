## @example simple_optimization/python/minimal_example.py
# This code contains a minimal example of an optimization problem that can be 
# built and solved using `panocpy`.

# %% Build the problem for PANOC+ALM (CasADi code, independent of panocpy)
import casadi as cs

name = "minimal_example"

# Make decision variables
x = cs.SX.sym("x")
y = cs.SX.sym("y")

# Make a parameter symbol
p = cs.SX.sym("p")

cost = (1 - x) ** 2 + p * (y - x ** 2) ** 2  # Rosenbrock function (parametrized by p)

constraint_g_cubic = (x - 1) ** 3 - y + 1
constraint_g_linear = x + y - 2

# Collect decision variables into one vector
X = cs.vertcat(x, y)

cost_function = cs.Function("f", [X, p], [cost])
g = cs.vertcat(constraint_g_cubic, constraint_g_linear)
g_function = cs.Function("g", [X, p], [g])

# %% Generate and compile C-code using `panocpy`
import panocpy as pa
import numpy as np

cgen, n, m, num_p = pa.generate_casadi_problem(name, cost_function, g_function)
# Code generator, dimension of decision variables, number of constraints (dual dimension), parameter dimension

# Compile and load the problem, and set the bounds
prob = pa.compile_and_load_problem(cgen, n, m, num_p, name)

prob.C.lowerbound = np.array([-1.5, -0.5])        # -1.5 <= x <= 1.5
prob.C.upperbound = np.array([1.5, 2.5])          # -0.5 <= y <= 2.5
prob.D.lowerbound = np.array([-np.inf, -np.inf])  # g_c <= 0
prob.D.upperbound = np.array([0, 0])              # g_l <= 0

# %% Construct a PANOC instance to serve as the inner solver

innersolver = pa.StructuredPANOCLBFGSSolver() # (with default parameters)

# %% Make an ALM solver with default parameters, using the PANOC solver

solver = pa.ALMSolver(pa.ALMParams(), innersolver)

# Set parameter to some value
prob.param = np.array([100.0])

# Set initial guesses at arbitrary values
x_sol = np.array([1.0, 2.0])
y_sol = np.zeros((m,))

# Solve the problem
x_sol, y_sol, stats = solver(prob, x_sol, y_sol)


print(stats["status"])

print(f"Obtained solution: {x_sol}")
print(f"Analytical solution: {(1., 1.)}")

# %% Plot the results

import matplotlib.pyplot as plt

x = np.linspace(-1.5, 1.5, 200)
y = np.linspace(-0.5, 2.5, 200)
X, Y = np.meshgrid(x, y)
Z = (1 - X) ** 2 + 100 * (Y - X ** 2) ** 2

plt.figure()

plt.contourf(X, Y, Z)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(1, 1, color="tab:red", label="Analytic")
plt.scatter(x_sol[0], x_sol[1], marker="x", color="tab:green", label="PANOC-ALM")
plt.legend()
plt.show()

# %%