from casadi import SX, Function, CodeGenerator, vertcat, jtimes, gradient

x = SX.sym("x")
y = SX.sym("y")
z = SX.sym("z")

v = SX.sym("v")
w = SX.sym("w")

# Formulate the NLP
f = x**2 + 100*z**2
g = z + (1-x)**2 - y

unknwns = vertcat(x, y, z)
objective = Function("f", [unknwns], [f], ["x"], ["f"])
grad_objective = Function("grad_f", [unknwns], [gradient(f, unknwns)], ["x"], ["grad_f"])
constraint = Function("g", [unknwns], [g], ["x"], ["g"])
grad_constraint = Function("grad_g_v", [unknwns, w], [jtimes(g, unknwns, w, True)], ["x", "w"], ["grad_g_v"])
cg = CodeGenerator("rosenbrock_functions.c")
cg.add(objective)
cg.add(grad_objective)
cg.add(constraint)
cg.add(grad_constraint)
cg.generate()
