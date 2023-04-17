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

# Define the bounds
C = [-0.25, -0.5], [1.5, 2.5]  # -0.25 <= x1 <= 1.5, -0.5 <= x2 <= 2.5
D = [-np.inf, -np.inf], [0, 0]  #         g1 <= 0,           g2 <= 0

# %% Generate and compile C-code for the objective and constraints using alpaqa
from alpaqa import minimize

problem = (
    minimize(f, x)  #       Objective function f(x)
    .subject_to_box(C)  #   Box constraints x ∊ C
    .subject_to(g, D)  #    General ALM constraints g(x) ∊ D
    .with_param(p, [1])  #  Parameter with default value (can be changed later)
).compile()

# You can change the bounds and parameters after loading the problem
problem.param = [10.0]
problem.D.lowerbound[1] = -1e20

# %% Build a solver with the default parameters
import alpaqa as pa

inner_solver = pa.PANOCSolver()
solver = pa.ALMSolver(inner_solver)

# %% Build a solver with custom parameters
inner_solver = pa.PANOCSolver(
    panoc_params={
        'max_iter': 1000,
        'stop_crit': pa.PANOCStopCrit.ApproxKKT,
        'print_interval': 1,
    },
    lbfgs_params={
        'memory': 10,
    },
)

solver = pa.ALMSolver(
    alm_params={
        'tolerance': 1e-10,
        'dual_tolerance': 1e-10,
        'initial_penalty': 50,
        'penalty_update_factor': 20,
        'print_interval': 1,
    },
    inner_solver=inner_solver
)

# %% Compute a solution

# Set initial guesses at arbitrary values
x0 = [0.1, 1.8]  # decision variables
y0 = [0.0, 0.0]  # Lagrange multipliers for g(x)

# Solve the problem
x_sol, y_sol, stats = solver(problem, x0, y0)

# Print the results
print(stats["status"])
print(f"Solution:      {x_sol}")
print(f"Multipliers:   {y_sol}")
print(f"Cost:          {problem.eval_f(x_sol):.5f}")
