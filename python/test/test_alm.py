def test_alm():
    import alpaqa as pa
    import alpaqa.casadi_loader as cl

    pp = pa.PANOCParams()
    lbfgs = pa.LBFGS({})
    panoc = pa.PANOCSolver(pp, lbfgs)
    solver = pa.ALMSolver(pa.ALMParams(), panoc)
    print(f"Solver: {solver} "
          f"({solver.__class__.__module__}.{solver.__class__.__qualname__})")

    import casadi as cs
    import numpy as np
    from pprint import pprint

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
    solver = pa.PANOCSolver(
        pa.PANOCParams(max_iter=200, print_interval=1),
        pa.LBFGS(pa.LBFGS.Params(memory=5)),
    )
    almparams = pa.ALMParams(max_iter=20, print_interval=1)
    almsolver = pa.ALMSolver(almparams, solver)
    pc = pa.with_counters(p)
    x0 = np.array([3, 3])
    y0 = np.zeros((m, ))
    x, y, stats = almsolver(pc, x=x0, y=y0)

    print("x", x)
    print("y", y)
    pprint(stats)
    print(pc.evaluations)
    assert stats['status'] == pa.SolverStatus.Converged

if __name__ == '__main__':
    test_alm()