import alpaqa as pa
import numpy as np
import casadi as cs

# Problem specification
# minimize  ½ xᵀHx
#  s.t.     Ax ≤ b
H = np.array([[3, -1], 
              [-1, 3]])
A = np.array([[2, 1]])
b = np.array([-1])

x = cs.SX.sym("x", 2)
f = cs.Function("f", [x], [0.5 * cs.dot(x, H @ x)])
g = cs.Function("g", [x], [A @ x])
problem = pa.generate_and_compile_casadi_problem(f, g, "example_name")

# Specify the bounds
problem.C.lowerbound = np.full((problem.n,), -np.inf)
problem.C.upperbound = np.full((problem.n,), np.inf)
problem.D.lowerbound = np.full((problem.m,), -np.inf)
problem.D.upperbound = b

# Settings for the outer augmented Lagrangian method
almparam = pa.ALMParams(
    ε = 1e-8,            # tolerance
    δ = 1e-8,
    Δ = 10,              # penalty update factor
    max_iter = 20,
    print_interval = 1,
)
# Settings for the inner PANOC solver
panocparam = pa.PANOCParams(
    max_iter = 500,
    print_interval = 10,
)
# Settings for the L-BFGS algorithm used by PANOC
lbfgsparam = pa.LBFGSParams(
    memory = 2,
)

# Create an ALM solver using PANOC as inner solver
solver = pa.ALMSolver(
    almparam,                               # Params for outer solver
    pa.PANOCSolver(panocparam, lbfgsparam), # Inner solver
)

# Initial guess
x = np.array([2, 2]) # decision variables
y = np.array([1])    # Lagrange multipliers

# Solve the problem
x_sol, y_sol, stats = solver(problem, x, y)

# Print the results
print(f"""
status: {stats['status']}
inner iterations: {stats['inner']['iterations']}
outer iterations: {stats['outer_iterations']}
elapsed time:     {stats['elapsed_time']}
x = {x_sol}
y = {y_sol}
f = {problem.f(x_sol)}
""")