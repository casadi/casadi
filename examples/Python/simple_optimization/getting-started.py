# %% Build the problem (CasADi code, independent of alpaqa)
import casadi as cs
import numpy as np

# Make symbolic decision variables
x1, x2 = cs.SX.sym("x1"), cs.SX.sym("x2")
# Collect decision variables into one vector
x = cs.vertcat(x1, x2)
# Make a parameter symbol
p = cs.SX.sym("p")

# Expressions for the objective function f and the constraints g
f_expr = (1 - x1) ** 2 + p * (x2 - x1 ** 2) ** 2
g_expr = cs.vertcat(
    (x1 - 0.5) ** 3 - x2 + 1,
    x1 + x2 - 1.5,
)

# Convert the symbolic expressions to CasADi functions
f = cs.Function("f", [x, p], [f_expr])
g = cs.Function("g", [x, p], [g_expr])

# Set the bounds
C = [-0.25, -0.5], [1.5, 2.5]  # -0.25 <= x1 <= 1.5, -0.5 <= x2 <= 2.5
D = [-np.inf, -np.inf], [0, 0]  #         g1 <= 0,           g2 <= 0
# Set the problem parameter
param = [10.0]

# %% Generate and compile C-code for the objective and constraints using alpaqa
import alpaqa.casadi_loader as cl

# Compile and load the problem
prob = cl.generate_and_compile_casadi_problem(
    f=f,  # minimize    f(x; param)
    C=C,  # subject to  x ∊ C
    g=g,  # subject to  g(x; param) ∊ D
    D=D,
    param=param,
)

# You can change the bounds and parameters after loading the problem
prob.param = [100.0]
prob.D.lowerbound[1] = -1e20

# %% Build a solver with the default parameters
import alpaqa as pa

innersolver = pa.PANOCSolver()
solver = pa.ALMSolver(innersolver)

# %% Build a solver with custom parameters
inner_solver = pa.PANOCSolver(
    panoc_params={
        "max_iter": 1000,
        "stop_crit": pa.PANOCStopCrit.ApproxKKT,
    },
    lbfgs_params={
        "memory": 10,
    },
)

solver = pa.ALMSolver(
    alm_params={
        "ε": 1e-10,
        "δ": 1e-10,
        "Σ_0": 0,
        "σ_0": 2,
        "Δ": 20,
        "print_interval": 1,
    },
    inner_solver=inner_solver,
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
print(f"Cost:          {prob.eval_f(x_sol):.5f}")
print(f"ε:             {stats['ε']}")
print(f"δ:             {stats['δ']}")
