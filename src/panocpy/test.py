import panocpy as pa
import numpy as np

solver = pa.PANOCSolver(pa.PANOCParams(), pa.LBFGSDirection(pa.LBFGSParams()))
assert str(solver) == "PANOCSolver<LBFGS>"
solver = pa.PANOCSolver(pa.PANOCParams(), pa.LBFGSParams())
assert str(solver) == "PANOCSolver<LBFGS>"


class Dir(pa.PANOCDirection):
    def __init__(self):
        super().__init__()

    def get_name(self):
        return "Dir"


assert str(Dir()) == "Dir"
solver = pa.PANOCSolver(pa.PANOCParams(), Dir())
assert str(solver) == "PANOCSolver<Dir>"

l = pa.LBFGSParams(cbfgs=pa.LBFGSParamsCBFGS(α=5))
assert l.cbfgs.α == 5
l.cbfgs.α = 100
assert l.cbfgs.α == 100

import casadi as cs

hess_prod = lambda L, x, v: cs.gradient(cs.jtimes(L, x, v, False), x)

n = 2
m = 2
x = cs.SX.sym("x", n)
λ = cs.SX.sym("λ", m)
v = cs.SX.sym("v", n)

Q = np.array([[1.5, 0.5], [0.5, 1.5]])
f_ = x.T @ Q @ x
g_ = x
L = f_ + cs.dot(λ, g_) if m > 0 else f_

f = cs.Function("f", [x], [f_])
grad_f = cs.Function("grad_f", [x], [cs.gradient(f_, x)])
g = cs.Function("g", [x], [g_])
grad_g_prod = cs.Function("grad_g_prod", [x, λ], [cs.jtimes(g_, x, λ, True)])
grad_gi = lambda x, i: grad_g_prod(x, np.eye(1, m, i))
Hess_L = cs.Function("Hess_L", [x, λ], [cs.hessian(L, x)[0]])
Hess_L_prod = cs.Function("Hess_L_prod", [x, λ, v], [hess_prod(L, x, v)])

p = pa.Problem(n, m)
p.f = f
p.grad_f = grad_f
p.g = g
p.grad_g_prod = grad_g_prod
p.grad_gi = grad_gi
p.hess_L = Hess_L
p.hess_L_prod = Hess_L_prod
p.D.lowerbound = [-np.inf, 0.5]
p.D.upperbound = [+np.inf, +np.inf]

x0 = np.array([3, 3])
y0 = np.zeros((m,))
Σ = 1e3 * np.ones((m,))
ε = 1e-8
solver = pa.PANOCSolver(
    pa.PANOCParams(max_iter=200, print_interval=1),
    pa.LBFGSDirection(pa.LBFGSParams(cbfgs=pa.LBFGSParamsCBFGS(α=2))),
)
x, y, err_z, stats = solver(p, Σ, ε, True, x0, y0)
print(x)
print(y)
print(err_z)
print(stats.ε)
print(stats.iterations)
print(stats.status)
print(stats.elapsed_time)
