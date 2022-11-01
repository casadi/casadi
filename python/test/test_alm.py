from copy import copy, deepcopy
import pickle
import alpaqa as pa
import numpy as np
from pprint import pprint

def test_alm():
    import alpaqa.casadi_loader as cl

    pp = pa.PANOCParams()
    lbfgs = pa.LBFGS({})
    panoc = pa.PANOCSolver(pp, lbfgs)
    solver = pa.ALMSolver(pa.ALMParams(), panoc)
    print(f"Solver: {solver} "
          f"({solver.__class__.__module__}.{solver.__class__.__qualname__})")

    import casadi as cs

    n = 2
    m = 2
    x = cs.SX.sym("x", n)

    Q = np.array([[1.5, 0.5], [0.5, 1.5]])
    f_ = 0.5 * x.T @ Q @ x
    g_ = x
    f = cs.Function("f", [x], [f_])
    g = cs.Function("g", [x], [g_])

    name = "testproblem"
    p = cl.generate_and_compile_casadi_problem(f, g, name=name)
    p = copy(p) # test copying/cloning
    p = deepcopy(p) # test copying/cloning
    print(p)
    print("C", p.C.lowerbound, p.C.upperbound)
    print("D", p.D.lowerbound, p.D.upperbound)
    p.D.lowerbound = [-np.inf, 0.5]
    p.D.upperbound = [+np.inf, +np.inf]
    solver = pa.PANOCSolver(
        pa.PANOCParams(max_iter=200, print_interval=1),
        pa.LBFGS(pa.LBFGS.Params(memory=5)),
    )
    almparams = pa.ALMParams(max_iter=20, print_interval=1, ε=1e-12, δ=1e-12)
    almsolver = pa.ALMSolver(almparams, solver)
    pc, counters = pa.problem_with_counters(p)
    x0 = np.array([3, 3])
    y0 = np.zeros((m, ))
    
    x, y, stats = almsolver(pc, x=x0, y=y0)

    print()
    print(counters)
    print(stats["status"])
    print("x", x)
    print("y", y)
    pprint(stats)
    assert stats['status'] == pa.SolverStatus.Converged

class MyProblem(pa.BoxConstrProblem):
    vec = np.ndarray
    Q = np.array([[1.5, 0.5], [0.5, 1.5]])
    A = np.eye(2)

    def __init__(self) -> None:
        super().__init__(self.Q.shape[0], self.A.shape[0])

    def eval_f(self, x: vec) -> float:
        return 0.5 * x.T @ self.Q @ x

    def eval_grad_f(self, x: vec, grad_fx: vec):
        grad_fx[:] = self.Q @ x

    def eval_g(self, x: vec, gx: vec):
        gx[:] = self.A @ x

    def eval_grad_g_prod(self, x: vec, y: vec, g: vec):
        g[:] = self.A.T @ y


def get_pickled_problem():
    p = MyProblem()
    p.D.lowerbound = [-np.inf, 0.5]
    p.D.upperbound = [+np.inf, +np.inf]
    return pickle.dumps(p)


def test_alm_inherit():

    p = pickle.loads(get_pickled_problem())
    solver = pa.PANOCSolver(
        pa.PANOCParams(max_iter=200, print_interval=1),
        pa.LBFGS(pa.LBFGS.Params(memory=5)),
    )
    almparams = pa.ALMParams(max_iter=20, print_interval=1, ε=1e-12, δ=1e-12)
    almsolver = pa.ALMSolver(almparams, solver)
    pc, counters = pa.problem_with_counters(p)
    x0 = np.array([3, 3])
    y0 = np.array([0, 0])
    x, y, stats = almsolver(pa.TEProblem(p), x=x0, y=y0)
    x, y, stats = almsolver(pc, x=x0, y=y0)

    print()
    print(counters)
    print(stats["status"])
    print("x", x)
    print("y", y)
    pprint(stats)
    assert stats['status'] == pa.SolverStatus.Converged


if __name__ == '__main__':
    test_alm()
    test_alm_inherit()
