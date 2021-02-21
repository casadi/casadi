from casadi import SX, Function, CodeGenerator, vertcat, jtimes, gradient
from sys import argv

if len(argv) < 2:
    print(f"Usage:    {argv[0]} <name>")
    exit(0)

x = SX.sym("x")
y = SX.sym("y")
z = SX.sym("z")
unknwns = vertcat(x, y, z)

w = SX.sym("w")

# Formulate the NLP
f = x**2 + 100*z**2
g = z + (1-x)**2 - y

cg = CodeGenerator(f"{argv[1]}.c")
cg.add(Function("f", [unknwns], 
                [f], 
                ["x"], ["f"]))
cg.add(Function("grad_f", [unknwns], 
                [gradient(f, unknwns)],
                ["x"], ["grad_f"]))
cg.add(Function("g", [unknwns],
                [g],
                ["x"], ["g"]))
cg.add(Function("grad_g", [unknwns, w],
                [jtimes(g, unknwns, w, True)],
                ["x", "w"], ["grad_g"]))
cg.generate()
