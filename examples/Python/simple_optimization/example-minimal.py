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
constr_function = cs.Function("g", [X, p], [g])

import numpy as np
import alpaqa
import alpaqa.casadi_loader as cl

# Generate and compile the optimized C code for the problem
prob = cl.generate_and_compile_casadi_problem(cost_function, constr_function)

prob.C.lowerbound = [-1.5, -0.5]        # -1.5 <= x <= 1.5
prob.C.upperbound = [1.5, 2.5]          # -0.5 <= y <= 2.5
prob.D.lowerbound = [-np.inf, -np.inf]  # g_c <= 0
prob.D.upperbound = [0, 0]              # g_l <= 0

# Construct a PANOC instance to serve as the inner solver
innersolver = alpaqa.PANOCSolver() # (with default parameters)
# Make an ALM solver with default parameters, using the PANOC solver
solver = alpaqa.ALMSolver(alpaqa.ALMParams(), innersolver)

# Solve the problem for the given parameter value
prob.param = [100.]
x_sol, y_sol, stats = solver(prob)

# Set initial guesses at arbitrary values
x_sol = np.array([1.0, 2.0])
y_sol = np.zeros((prob.m,))

# Solve the problem
x_sol, y_sol, stats = solver(prob, x_sol, y_sol)

print(stats["status"])
print(f"Solution: {x_sol}")
print(f"ε:        {stats['ε']}")
print(f"δ:        {stats['δ']}")
