# %% alpaqa nonlinear regression example

import alpaqa as pa
import alpaqa.casadi_loader as cl
import casadi as cs
import numpy as np
from pprint import pprint

# %% Build the problem (CasADi code, independent of alpaqa)

# Symbolic model parameters
p = cs.MX.sym("params", 3)
# Symbolic data vectors
N_data = 50
data_x = cs.MX.sym("data_x", N_data)
data_y = cs.MX.sym("data_y", N_data)
data = cs.horzcat(data_x, data_y)

# Objective function is the squared error between the model and the data
x, y = cs.MX.sym("x"), cs.MX.sym("y")
model = cs.Function("model", [x, p], [p[0] * cs.sin(p[1] * x) + p[2]])
sample_error = cs.Function("err", [x, y, p], [y - model(x, p)])
sum_sq_error = cs.sumsqr(sample_error.map(N_data)(data_x.T, data_y.T, p))

# Convert the symbolic expression to a CasADi function
f = cs.Function("f", [p, cs.vec(data)], [sum_sq_error])

# %% Generate and compile C-code for the objective and constraints using alpaqa

# Compile and load the problem (without general constraints)
prob = cl.generate_and_compile_casadi_problem(f, None, second_order="psi_prod")
# Optionally, add constraints on the parameters
prob.C.lowerbound[1] = 3
prob.C.upperbound[1] = 7

# %% Generate some data

true_params = [6, 5, 4]
true_model = np.vectorize(lambda x: model(x, true_params))
rng = np.random.default_rng(12345)
data_x = np.linspace(-1, +1, N_data, endpoint=True)
data_y = true_model(data_x) + rng.standard_normal(N_data)

# Add data to the problem
prob.param = np.concatenate((data_x, data_y))

# %% Solve the problem using alpaqa's PANTR solver

solver = pa.PANTRSolver({"print_interval": 1})
# Add evaluation counters to the problem
cnt = pa.problem_with_counters(prob)
sol_params, stats = solver(cnt.problem, {"tolerance": 1e-10})

# %% Print the results

print(f"\nSolution:      {sol_params}")
print(stats["status"])
pprint(stats)
print("\nEvaluations:")
print(cnt.evaluations)

# %% Plot the results

import matplotlib.pyplot as plt

model_str = r"${:.2f}\ \sin({:.2f} x) + {:.2f}$"
sol_model = np.vectorize(lambda x: model(x, sol_params))
plt.figure(figsize=(8, 5))
x_fine = np.linspace(-1, 1, 256)
plt.plot(data_x, data_y, "x", label="Samples")
plt.plot(x_fine, true_model(x_fine),
         label="True model\n" + model_str.format(*true_params))
plt.plot(x_fine, sol_model(x_fine), "--",
         label="Est. model\n" + model_str.format(*sol_params))
plt.legend(loc='upper right')
plt.title("PANTR nonlinear regression example")
plt.tight_layout()
plt.show()

# %%
