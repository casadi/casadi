from typing import assert_type

import casadi as ca

x = ca.MX.sym("x")
y = ca.sin(x)
assert_type(x, ca.MX)
assert_type(y, ca.MX)
assert_type(ca.vertcat(x, y), ca.MX)

sx = ca.SX.sym("sx")
assert_type(ca.gradient(ca.dot(sx, sx), sx), ca.SX)

dense = ca.DM.eye(3)
assert_type(dense, ca.DM)
assert_type(ca.mtimes(dense, dense), ca.DM)

f = ca.Function("f", [x], [y])
assert_type(f, ca.Function)

solver = ca.nlpsol("solver", "ipopt", {"x": x, "f": y}, {})
assert_type(solver, ca.Function)

opti = ca.Opti()
decision = opti.variable()
assert_type(decision, ca.MX)
assert_type(opti.solve(), ca.OptiSol)
