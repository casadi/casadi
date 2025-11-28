# This simple MIP comes from https://docs.gurobi.com/projects/examples/en/current/examples/python/mip1.html#subsubsectionmip1-py
# Here it is re-written to match the Casadi syntax

import casadi as ca

# Binary decision variables
x = ca.SX.sym("x")
y = ca.SX.sym("y")
z = ca.SX.sym("z")

# Constraints
g = []
g.append(x + 2*y + 3*z)
lbg = [-ca.inf]
ubg = [4]
g.append(x + y)
lbg.append(1)
ubg.append(ca.inf)

# Objective: minimize -(x + y + 2 z)
f = -(x + y + 2*z)

# Standard behavior of gurobi -- printing to console is allowed
# The option suppress_all_output prevents gurobi to print to console.
print("\n\n#################### STANDARD GUROBI BEHAVIOR ####################")
solver = ca.qpsol('solver', 'gurobi',
                  {'f': f, 'g': ca.vertcat(*g), 'x': ca.vertcat(x, y, z)},
                  {'discrete': [1, 1, 1]})

# Solve
sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")

sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")

# Even if you set gurobi options OutputFlag or LogToConsole some information like the license banner will still be printed, making the log cluttered in case of repetitive optimization like in model predictive control
print("\n\n#################### GUROBI WITH OutputFlag=0 ####################")
solver = ca.qpsol('solver', 'gurobi',
                  {'f': f, 'g': ca.vertcat(*g), 'x': ca.vertcat(x, y, z)},
                  {'discrete': [1, 1, 1],
                   "gurobi.OutputFlag": 0})

# Solve
sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")

sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")

# The new option suppress_all_output prevents gurobi to print to console.
# It follows the guidelines of setting Gurobi OutputFlag before creating the gurobi environment like reported here:
# https://support.gurobi.com/hc/en-us/articles/360044784552-How-do-I-suppress-all-console-output-from-Gurobi
print("\n\n#################### GUROBI WITH new flag suppress_all_output=True ####################")
solver = ca.qpsol('solver', 'gurobi',
                  {'f': f, 'g': ca.vertcat(*g), 'x': ca.vertcat(x, y, z)},
                  {'discrete': [1, 1, 1], 'suppress_all_output': True})

# Solve
sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")

sol = solver(lbx=0, ubx=1, lbg=lbg, ubg=ubg)
print(f"Optimal solution: {sol['x'].full().squeeze()}")
