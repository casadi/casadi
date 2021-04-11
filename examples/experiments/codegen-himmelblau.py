from casadi import SX, Function, CodeGenerator, vertcat, jtimes, gradient, hessian
from sys import argv

if len(argv) < 2:
    print(f"Usage:    {argv[0]} <name>")
    exit(0)

x = SX.sym("x")
y = SX.sym("y")
unknwns = vertcat(x, y)

w = SX.sym("w", 2)
v = SX.sym("v", 2)
λ = SX.sym("λ", 2)

# Formulate the NLP
f = (x**2 + y - 11)**2 + (x + y**2 - 7)**2
g1 = (x-3)**2 + (y-2)**2
g2 = x**2
g = vertcat(g1, g2)
L = f + λ.T @ g

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
cg.add(Function("hess_L", [unknwns, λ],
                [hessian(L, unknwns)[0]],
                ["x", "y"], ["hess_L"]))
cg.add(Function("hess_L_prod", [unknwns, λ, v],
                [gradient(jtimes(L, unknwns, v, False), unknwns)],
                ["x", "y", "v"], ["hess_L_prod"]))
cg.generate()