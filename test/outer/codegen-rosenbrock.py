from sys import argv
import casadi as cs
import alpaqa.casadi_generator as cg

if len(argv) < 2:
    print(f"Usage:    {argv[0]} <name>")
    exit(0)

x = cs.SX.sym("x")
y = cs.SX.sym("y")
unknowns = cs.vertcat(x, y)

p = cs.SX.sym("p")

# Formulate the NLP
# (https://web.archive.org/web/20201111225115/https://web.casadi.org/blog/opti/)
f = (1 - x)**2 + (y - x**2)**2
g = x**2 + p * y**2

C = [0, 0], [5, 5]
D = [0], [1]
param = [2]

cgen = cg.generate_casadi_problem(
    cs.Function("f", [unknowns, p], [f]),
    cs.Function("g", [unknowns, p], [g]),
    second_order="full",
    name=argv[1],
)
cgen.generate()
cg.write_casadi_problem_data(argv[1], C, D, param)
