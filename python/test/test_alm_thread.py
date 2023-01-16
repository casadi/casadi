from copy import deepcopy
import alpaqa as pa
import numpy as np
import concurrent.futures
import pytest
import os


def test_alm_threaded():
    import alpaqa.casadi_loader as cl

    pp = pa.PANOCParams(max_no_progress=100, max_iter=100)
    lbfgs = pa.LBFGSDirection()
    panoc = pa.PANOCSolver(pp, lbfgs)
    alm_params = pa.ALMParams(ε=1e-200, δ=1e-200, max_iter=200, print_interval=0)
    solver = pa.ALMSolver(alm_params, panoc)

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
    p.D.lowerbound = [-np.inf, 0.5]
    p.D.upperbound = [+np.inf, +np.inf]
    p = pa.Problem(p)

    def good_experiment():
        _, _, stats = deepcopy(solver)(deepcopy(p), asynchronous=True)
        return stats

    def bad_experiment1():
        _, _, stats = solver(deepcopy(p), asynchronous=True)
        return stats

    def bad_experiment2():
        _, _, stats = deepcopy(solver)(p, asynchronous=True)
        return stats

    def run(experiment):
        with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
            futures = (pool.submit(experiment) for _ in range(200))
            for future in concurrent.futures.as_completed(futures):
                stats = future.result()
                assert stats["status"] == pa.SolverStatus.MaxIter

    run(good_experiment)
    with pytest.raises(
        RuntimeError, match=r"^Same instance of ALMSolver<PANOCSolver<LBFGS"
    ) as e:
        run(bad_experiment1)
    print(e.value)
    with pytest.raises(
        RuntimeError, match=r"^Same instance of type alpaqa::TypeErasedProblem"
    ) as e:
        run(bad_experiment2)
    print(e.value)


if __name__ == "__main__":
    test_alm_threaded()
    print("done.")
