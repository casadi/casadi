# %% alpaqa lasso example

import alpaqa as pa
from jax.config import config
import jax.numpy as jnp
from jax import grad, jit
from jax import random
from pprint import pprint

config.update("jax_enable_x64", True)

scale = 5000
n, m = scale, scale * 2
sparsity = 0.02
key = random.PRNGKey(0)

# %% Generate some data

key, *subkeys = random.split(key, 5)
# Random data matrix A
A = random.uniform(subkeys[0], (m, n), minval=-1, maxval=1)
# Sparse solution x_exact
x_exact = random.uniform(subkeys[1], (n,), minval=-0.1, maxval=1)
x_exact_zeros = random.uniform(subkeys[2], (n,), minval=0, maxval=1) > sparsity
x_exact = x_exact.at[x_exact_zeros].set(0)
# Noisy right-hand side b
b = A @ x_exact + 0.1 * random.normal(subkeys[3], (m,))

# %% Build the problem

# Quadratic loss plus l1-regularization
# minimize  ½‖Ax - b‖² + λ‖x‖₁

λ = 0.0025 * m


def loss(x):
    err = A @ x - b
    return 0.5 * jnp.dot(err, err)


class LassoProblem(pa.BoxConstrProblem):
    def __init__(self):
        super().__init__(n, 0)
        self.C.lowerbound[:] = 0  # Positive lasso
        self.l1_reg = [λ]  # Regularization
        self.jit_loss = jit(loss)
        self.jit_grad_loss = jit(grad(loss))

    def eval_f(self, x):  # Cost function
        return self.jit_loss(x)

    def eval_grad_f(self, x, grad_f):  # Gradient of the cost
        grad_f[:] = self.jit_grad_loss(x)


prob = LassoProblem()

# %% Solve the problem using alpaqa's PANOC solver

opts = {
    "max_iter": 100,
    "stop_crit": pa.FPRNorm,
    # Use a laxer tolerance because large problems have more numerical errors:
    "quadratic_upperbound_tolerance_factor": 1e-12,
}
direction = pa.StructuredLBFGSDirection({"memory": 5}, {"hessian_vec": False})
# direction = pa.LBFGSDirection({"memory": 5})
solver = pa.PANOCSolver({"print_interval": 5} | opts, direction)
# Add evaluation counters to the problem
cnt = pa.problem_with_counters(prob)
# Solve the problem
sol, stats = solver(cnt.problem, {"tolerance": 1e-10})

# %% Print the results

final_f = prob.eval_f(sol)
print()
pprint(stats)
print()
print("Evaluations:")
print(cnt.evaluations)
print(f"Cost:          {final_f + stats['final_h']}")
print(f"Loss:          {final_f}")
print(f"Regularizer:   {stats['final_h']}")
print(f"FP Residual:   {stats['ε']}")
print(f"Run time:      {stats['elapsed_time']}")
print(stats["status"])

# %% Plot the results

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 5))
plt.plot(x_exact, ".-", label="True solution")
plt.plot(sol, ".-", label="Estimated solution")
plt.legend()
plt.title("PANOC lasso example: solution")
plt.tight_layout()
plt.figure(figsize=(8, 5))
plt.plot(A @ x_exact, ".-", label="True solution")
plt.plot(A @ sol, ".-", label="Estimated solution")
plt.legend()
plt.title("PANOC lasso example: right-hand side")
plt.tight_layout()
plt.show()
