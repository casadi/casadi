#%%

print("starting tests")

import alpaqa as pa
import numpy as np
from pprint import pprint

print("1")
solver = pa.PANOCSolver(pa.PANOCParams(), pa.LBFGSDirection(pa.LBFGSParams()))
print("2")
assert str(solver) == "PANOCSolver<LBFGS>"
print("3")
solver = pa.PANOCSolver(pa.PANOCParams(), pa.LBFGSParams())
print("4")
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
    pa.LBFGSParams(memory=5),
)
x, y, err_z, stats = solver(p, Σ, ε, x0, y0)
print(x)
print(y)
print(err_z)
pprint(stats)

solver = pa.PANOCSolver(
    pa.PANOCParams(max_iter=200, print_interval=1),
    pa.LBFGSParams(memory=5),
)
almparams = pa.ALMParams(max_iter=20, print_interval=1, preconditioning=False)
almsolver = pa.ALMSolver(almparams, solver)
x, y, stats = almsolver(p, x=x0, y=y0)

print(x)
print(y)
pprint(stats)

solver = pa.StructuredPANOCLBFGSSolver(
    pa.StructuredPANOCLBFGSParams(max_iter=200, print_interval=1),
    pa.LBFGSParams(memory=5),
)
almparams = pa.ALMParams(max_iter=20, print_interval=1, preconditioning=False)
almsolver = pa.ALMSolver(almparams, solver)
x, y, stats = almsolver(p, x=x0, y=y0)


class CustomInnerSolver(pa.InnerSolver):
    def __init__(self):
        super().__init__()
        self.solver = pa.PANOCSolver(
            pa.PANOCParams(max_iter=200, print_interval=1),
            pa.LBFGSParams(memory=5),
        )

    def get_name(self):
        return self.solver.get_name()

    def stop(self):
        return self.solver.stop()

    def __call__(self, problem, Σ, ε, always_overwrite_results, x, y):
        # TODO: always_overwrite_results
        x, y, err_z, stats = self.solver(problem, Σ, ε, x, y)

        def accumulate(acc: dict, s: dict):
            for k, v in s.items():
                if not k in ["status", "ε", "accumulator"]:
                    acc[k] = acc[k] + v if k in acc else v

        stats["accumulator"] = {"accumulate": accumulate}
        return x, y, err_z, stats


solver = CustomInnerSolver()
almparams = pa.ALMParams(max_iter=20, print_interval=1, preconditioning=False)
almsolver = pa.ALMSolver(almparams, solver)
x, y, stats = almsolver(p, x=x0, y=y0)

print(x)
print(y)
pprint(stats)

try:
    old_x0 = x0
    x0 = np.zeros((666,))
    sol = almsolver(p, x=x0, y=y0)
except ValueError as e:
    assert e.args[0] == "Length of x does not match problem size problem.n"

x0 = old_x0

# %%

n = 2
m = 2
x = cs.SX.sym("x", n)
λ = cs.SX.sym("λ", m)
v = cs.SX.sym("v", n)

Q = np.array([[1.5, 0.5], [0.5, 1.5]])
f_ = 0.5 * x.T @ Q @ x
g_ = x
f = cs.Function("f", [x], [f_])
g = cs.Function("g", [x], [g_])

name = "testproblem"
p = pa.generate_and_compile_casadi_problem(f, g, name=name)
p.D.lowerbound = [-np.inf, 0.5]
p.D.upperbound = [+np.inf, +np.inf]
solver = pa.StructuredPANOCLBFGSSolver(
    pa.StructuredPANOCLBFGSParams(max_iter=200, print_interval=1),
    pa.LBFGSParams(memory=5),
)
almparams = pa.ALMParams(max_iter=20, print_interval=1, preconditioning=False)
almsolver = pa.ALMSolver(almparams, solver)
x, y, stats = almsolver(p, x=x0, y=y0)

print(x)
print(y)
pprint(stats)

# %%

n = 2
m = 2
x = cs.SX.sym("x", n)
p = cs.SX.sym("p", 3)

p0 = np.array([1.5, 0.5, 1.5])

Q = cs.vertcat(cs.horzcat(p[0], p[1]), cs.horzcat(p[1], p[2]))
f_ = 0.5 * x.T @ Q @ x
g_ = x
f = cs.Function("f", [x, p], [f_])
g = cs.Function("g", [x, p], [g_])

name = "testproblem"
prob = pa.generate_and_compile_casadi_problem(f, g, name=name)
prob.D.lowerbound = [-np.inf, 0.5]
prob.D.upperbound = [+np.inf, +np.inf]
prob.param = p0
solver = pa.StructuredPANOCLBFGSSolver(
    pa.StructuredPANOCLBFGSParams(max_iter=200, print_interval=1),
    pa.LBFGSParams(memory=5),
)
almparams = pa.ALMParams(max_iter=20, print_interval=1, preconditioning=False)
almsolver = pa.ALMSolver(almparams, solver)
x, y, stats = almsolver(prob, x=x0, y=y0)

print(x)
print(y)
pprint(stats)

# %% Make sure that the problem is copied

prob.param = [1, 2, 3]
assert np.all(prob.param == [1, 2, 3])
prob1 = pa.ProblemWithParamWithCounters(prob)
print(prob.param)
print(prob1.param)
assert np.all(prob.param == [1, 2, 3])
assert np.all(prob1.param == [1, 2, 3])
prob1.param = [42, 43, 44]
print(prob.param)
print(prob1.param)
assert np.all(prob.param == [1, 2, 3])
assert np.all(prob1.param == [42, 43, 44])
print(prob.f([1, 2]))
print(prob1.f([1, 2]))
assert prob.f([1, 2]) == 21 / 2
assert prob1.f([1, 2]) == 390 / 2
assert prob1.evaluations.f == 2

prob2 = pa.ProblemWithCounters(prob) # params are not copied!
print(prob.f([1, 2]))
print(prob2.f([1, 2]))
assert prob.f([1, 2]) == 21 / 2
assert prob2.f([1, 2]) == 21 / 2
prob.param = [2, 1, 3]
print(prob.f([1, 2]))
print(prob2.f([1, 2]))
assert prob.f([1, 2]) == 18 / 2
assert prob2.f([1, 2]) == 18 / 2
assert prob1.evaluations.f == 2
assert prob2.evaluations.f == 4

# %%

prob1.param = p0
x, y, stats = almsolver(prob1)  # without initial guess

print(x)
print(y)
pprint(stats)
print(prob1.evaluations.f)
print(prob1.evaluations.grad_f)
print(prob1.evaluations.g)
print(prob1.evaluations.grad_g_prod)

# %%

f = lambda x: float(np.cosh(x) - x * x + x)
grad_f = lambda x: np.sinh(x) - 2 * x + 1
C = pa.Box([10], [-2.5])
x0 = [5]
x, stats = pa.panoc(
    f, grad_f, C, x0, 1e-12, pa.PANOCParams(print_interval=1), pa.LBFGSParams()
)
print(x)
pprint(stats)

# %%

f = lambda x: float(np.cosh(x) - x * x + x)
grad_f = lambda x: np.sinh(x) - 2 * x + 1
C = pa.Box([10], [-2.5])
x, stats = pa.panoc(f, grad_f, C, params=pa.PANOCParams(print_interval=1))
print(x)
pprint(stats)

# %%

try:
    pa.PANOCParams(max_iter=1e3)
    assert False
except RuntimeError as e:
    print(e)

# %%

print("Success!")