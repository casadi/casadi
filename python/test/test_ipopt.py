import alpaqa as pa
import numpy as np
from pprint import pprint
import pytest


@pytest.mark.skipif(
    not pa.with_casadi or not pa.with_ipopt, reason="requires CasADi and Ipopt"
)
def test_ipopt_pyapi_compile():
    import casadi as cs

    n = 2
    m = 2
    x = cs.SX.sym("x", n)

    Q = np.array([[1.5, 0.5], [0.5, 1.5]])
    f = 0.5 * x.T @ Q @ x
    g = x
    D = [-np.inf, 0.5], [+np.inf, +np.inf]
    p = pa.minimize(f, x).subject_to(g, D).compile(second_order="L")
    print(p)
    solver = pa.IpoptSolver({})
    cnt = pa.problem_with_counters(p)
    x0 = np.array([3, 3])
    y0 = np.zeros((m,))

    x, y, stats = solver(cnt.problem, x=x0, y=y0)

    print()
    print(cnt.evaluations)
    print(stats["status"])
    print("x", x)
    print("y", y)
    pprint(stats)
    assert stats["status"] == pa.IpoptSolverReturn.SUCCESS
    assert np.linalg.norm(x - [-1 / 6, 0.5]) < 1e-5
    assert np.linalg.norm(y - [0, -2 / 3]) < 1e-5
