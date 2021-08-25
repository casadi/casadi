## @example simple_optimization/python/minimal_example.py
# This code contains a minimal example of an optimization problem that can be built and solved using `panocpy`. 
import casadi as cs 
import panocpy as pa  
def compile_and_load_problem(cgen: cs.CodeGenerator, n: int, m: int, name: str = "PANOC_ALM_problem", ) -> pa.Problem:
    """Compile the C-code using the given code-generator and load it as a 
    panocpy Problem. 

    Args:
        cgen (cs.CodeGenerator): Code generator to generate C-code for the costs and the constraints with. 
        n (int): Dimensions of the decision variables (primal dimension) 
        m (int): Number of nonlinear constraints (dual dimension)
        name (str, optional): String description of the problem. Defaults to "PANOC_ALM_problem".
    """
    
    with TemporaryDirectory(prefix="") as tmpdir:
        cfile = cgen.generate(tmpdir)
        sofile = os.path.join(tmpdir, f"{name}.so")
        os.system(f"cc -fPIC -shared -O3 -march=native {cfile} -o {sofile}")
        prob = pa.load_casadi_problem_with_param(sofile, n, m)
    return prob 


#%% Import necessary libraries and generate the MPC problem

import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import timedelta
import os
import sys

sys.path.append(os.path.dirname(__file__))

# %% Build the problem for PANOC+ALM (CasADi code, independent of panocpy)

name = "minimal_example"

# Make decision variables 
x = cs.SX.sym("x")
y = cs.SX.sym("y")

# Make a parameter symbol 
p = cs.SX.sym("p")

cost = (1 - x)**2 + p * (y -x**2)**2  # Rosenbrock function (parametrized by p) 

constraint_g_cubic = (x-1)**3 - y + 1
constraint_g_linear = x + y - 2

# Collect decision variables into one vector 
X = cs.vertcat(x,y)

cost_function = cs.Function("f", [X, p], [cost])
g = cs.vertcat(constraint_g_cubic, constraint_g_linear)
g_function = cs.Function("g", [X, p], [g])

ipopt_solver = cs.nlpsol("solver", "ipopt", {"f": cost, "g": g, "x": X, "p": p})
solution_ipopt = ipopt_solver(p=100.)
x_ipopt = solution_ipopt["x"]


# %% Generate and compile C-code using `panocpy`
import panocpy as pa
from tempfile import TemporaryDirectory

cgen, n, m, num_p = pa.generate_casadi_problem(name, cost_function, g_function)
# Code generator, dimension of decision variables, number of constraints (dual dimension), parameter dimension 

# Compile and load the problem, and set the bounds 
prob = compile_and_load_problem(cgen, n, m, name)

prob.C.lowerbound = np.array([-1.5, -0.5])  # -1.5 <= x <= 1.5 
prob.C.upperbound = np.array([ 1.5,  2.5])  # -0.5 <= y <= 2.5
prob.D.lowerbound = np.array([-np.inf, -np.inf])  # g_c <= 0 
prob.D.upperbound = np.array([0, 0])              # g_l <= 0

#%% Construct a PANOC instance to serve as the inner solver (with default parameters) 

innersolver = pa.PANOCSolver(pa.PANOCParams(), pa.LBFGSParams())

#%% Make an ALM solver with default parameters 

almparams = pa.ALMParams()
solver = pa.ALMSolver(almparams, innersolver)

# Set parameter to some value 
prob.param = np.array([100.])

# set initial guesses at arbitrary values 
x_sol = np.array([1., 2.]) 
y_sol = np.zeros((m,))

y_sol, x_sol, stats = solver(prob, y_sol, x_sol)


print(stats["status"])

print(f"Obtained solution: {x_sol}")
# print(f"IPOPT solution: {x_ipopt}")
print(f"Analytical solution: {(1., 1.)}")


import matplotlib.pyplot as plt 

x = np.linspace(-1.5, 1.5, 200)
y = np.linspace(-0.5, 2.5, 200)
X,Y = np.meshgrid(x,y)
Z = (1 - X)**2 + 100 * (Y - X**2)**2 

plt.figure() 

plt.contourf(X,Y,Z) 
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(1, 1, color="tab:red", label="Analytic")
plt.scatter(x_sol[0], x_sol[1], marker="x", color="tab:green", label="PANOC-ALM")
# x_ipopt = x_ipopt.toarray()
# plt.scatter(x_ipopt[0,0], x_ipopt[1,0], marker="x", color="tab:blue", label="IPOPT")
plt.legend()
plt.show()
# %%
