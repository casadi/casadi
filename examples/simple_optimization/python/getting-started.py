## @example simple_optimization/python/getting-started.py
# This code contains a minimal example of an optimization problem that can be 
# built and solved using `alpaqa`.

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

# %% Build a solver with the default parameters
innersolver = pa.StructuredPANOCLBFGSSolver()
solver = pa.ALMSolver(innersolver)

# %% Build a solver with custom parameters
inner_solver = pa.StructuredPANOCLBFGSSolver(
    panoc_params={
        'max_iter': 1000,
        'stop_crit': pa.PANOCStopCrit.ApproxKKT,
    },
    lbfgs_params={
        'memory': 10,
    },
)

solver = pa.ALMSolver(
    alm_params={
        'ε': 1e-10,
        'δ': 1e-10,
        'Σ_0': 0,
        'σ_0': 2,
        'Δ': 20,
    },
    inner_solver=inner_solver
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
print(f"ε:             {stats['ε']}")
print(f"δ:             {stats['δ']}")
