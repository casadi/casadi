# %% alpaqa lasso example

import alpaqa as pa
import alpaqa.casadi_loader as cl
import casadi as cs
import numpy as np
from pprint import pprint

scale = 50
n, m = scale, scale * 2
sparsity = 0.2
rng = np.random.default_rng(0)

# %% Build the problem (CasADi code, independent of alpaqa)

# Quadratic loss plus l1-regularization
# minimize  ‖Ax - b‖² + λ‖x‖₁

A = rng.uniform(-1, 1, (m, n))
x_exact = rng.uniform(0, 1, n)
x_exact[rng.uniform(0, 1, n) > sparsity] = 0
b = A @ x_exact + rng.normal(0, 0.1, m)
λ = 5

# Symbolic solution
x = cs.MX.sym("x", n)
# Objective function is squared norm of Ax - b
f = cs.Function("f", [x], [cs.sumsqr(A @ x - b)])

# %% Generate and compile C-code for the objective and constraints using alpaqa

# Compile and load the problem (without general constraints)
prob = cl.generate_and_compile_casadi_problem(f, None, sym=cs.MX.sym)
# Set the regularization parameter
prob.l1_reg = [λ]

# %% Solve the problem using alpaqa's PANOC solver

direction = pa.LBFGSDirection({"memory": scale})
solver = pa.PANOCSolver({"print_interval": 10}, direction)
# Add evaluation counters to the problem
cnt = pa.problem_with_counters(prob)
sol, stats = solver(cnt.problem, {"tolerance": 1e-10})

# %% Print the results

print(f"\nSolution:      {sol}")
print(stats["status"])
pprint(stats)
print("\nEvaluations:")
print(cnt.evaluations)

# %% Plot the results

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 5))
plt.plot(x_exact, ".-", label="True solution")
plt.plot(sol, ".-", label="Estimated solution")
plt.legend()
plt.title("PANOC lasso example")
plt.tight_layout()
plt.show()
