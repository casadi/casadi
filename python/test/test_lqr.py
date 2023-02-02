from dataclasses import dataclass
import alpaqa as pa
import casadi as cs
import numpy as np
import numpy.random as nprand
import numpy.linalg as la
import alpaqa.casadi_loader as cl
import scipy.sparse as sp
import os

@dataclass
class OCProblemData:
    f: cs.Function
    h: cs.Function
    h_N: cs.Function
    l: cs.Function
    l_N: cs.Function
    c: cs.Function
    N: int

def convert_to_single_shooting(ocp: OCProblemData):
    nx = ocp.f.size1_in(0)
    nu = ocp.f.size1_in(1)
    mpc_x0 = cs.SX.sym("x0", nx)
    mpc_u = cs.SX.sym("u", (nu, ocp.N))  # Inputs
    mpc_x = ocp.f.mapaccum(ocp.N)(mpc_x0, mpc_u)  # Simulated states

    mpc_h = ocp.h.map(ocp.N)(cs.horzcat(mpc_x0, mpc_x[:, :-1]), mpc_u)
    mpc_stage_costs = cs.sum2(ocp.l.map(ocp.N)(mpc_h))
    mpc_terminal_cost = ocp.l_N(ocp.h_N(mpc_x[:, -1]))
    mpc_cost = cs.Function('f_mpc', [cs.vec(mpc_u), mpc_x0],
        [mpc_stage_costs + mpc_terminal_cost])

    g = [ocp.c(mpc_x0)]
    g += [ocp.c(mpc_x[:, i]) for i in range(ocp.N - 1)]
    # TODO: separate terminal constraint
    g += [ocp.c(mpc_x[:, ocp.N - 1])]
    mpc_constr = cs.Function("g", [cs.vec(mpc_u), mpc_x0],
        [cs.vertcat(*g)])

    return mpc_cost, mpc_constr

def convert_to_multiple_shooting_alm(ocp: OCProblemData, y, μ, D, D_N):
    N = ocp.N
    nx = ocp.f.size1_in(0)
    nu = ocp.f.size1_in(1)
    mpc_u = cs.SX.sym("u", (nu, N))  # Inputs
    mpc_x = cs.SX.sym("x", (nx, N + 1)) # States

    mpc_h = ocp.h.map(N)(mpc_x[:, :-1], mpc_u)
    mpc_stage_costs = cs.sum2(ocp.l.map(N)(mpc_h))
    mpc_terminal_cost = ocp.l_N(ocp.h_N(mpc_x[:, -1]))
    mpc_xu = cs.vertcat(cs.vec(cs.vertcat(mpc_x[:, :-1], mpc_u)), mpc_x[:, -1])

    penalty = 0
    for i in range(N):
        proj_D = lambda g: cs.fmax(D.lowerbound, cs.fmin(g, D.upperbound))
        ζ = ocp.c(mpc_x[:, i]) + y[i * nx:i * nx + nx] / μ[i]
        penalty += 0.5 * μ[i] * cs.sum1((ζ - proj_D(ζ))**2)
    proj_D_N = lambda g: cs.fmax(D_N.lowerbound, cs.fmin(g, D_N.upperbound))
    ζ = ocp.c(mpc_x[:, N]) + y[N * nx:N * nx + nx] / μ[N]
    penalty += 0.5 * μ[N] * cs.sum1((ζ - proj_D_N(ζ))**2)

    mpc_cost = cs.Function('f_mpc', [mpc_xu],
        [mpc_stage_costs + mpc_terminal_cost + penalty])

    return mpc_cost

def test_lqr():
    if not pa.with_casadi_ocp:
        return
    tol = 1e-10
    rng = nprand.default_rng(112233)

    # States and inputs
    N = 51
    nx, nu = 5, 7
    nc = nc_N = nx
    x_sym = cs.SX.sym("x", nx)
    u_sym = cs.SX.sym("u", nu)

    # Dynamics
    x_next = cs.vertcat(
        cs.tanh(x_sym[0]) + u_sym[0] * u_sym[1],
        cs.tanh(x_sym[1]) * cs.tanh(x_sym[0]) + u_sym[1] * u_sym[2],
        cs.tanh(x_sym[2]) * cs.tanh(x_sym[1]) + u_sym[3] * u_sym[4],
        cs.tanh(x_sym[3]) * cs.tanh(x_sym[2]) + u_sym[1] * u_sym[5],
        cs.tanh(x_sym[4]) * cs.tanh(x_sym[3]) + u_sym[2] * u_sym[6],
    )
    f = cs.Function("f", [x_sym, u_sym],
        [x_next])
    # Output
    h = cs.Function("h", [x_sym, u_sym],
        [cs.vertcat(cs.sinh(x_sym), 1. / (u_sym**2 + 1.))])
    hN = cs.Function("hN", [x_sym],
        [1. / (x_sym**2 + 1.)])
    # Cost
    L = rng.uniform(-1, 1, (nx + nu, nx + nu))
    L = L.T @ L
    hsym = cs.SX.sym("h", nx + nu)
    l = cs.Function("l", [hsym], 
        [cs.log(cs.cosh(hsym).T @ L @ cs.cosh(hsym))])
    LN = rng.uniform(-1, 1, (nx, nx))
    LN = LN.T @ LN
    hNsym = cs.SX.sym("hN", nx)
    lN = cs.Function("lN", [hNsym],
        [cs.log(cs.cosh(hNsym).T @ LN @ cs.cosh(hNsym))])
    # Constraints
    c = cs.Function("c", [x_sym],
        [cs.sinh(x_sym)])

    # Compile problem
    problem = cl.generate_and_compile_casadi_control_problem(
        N, f=f, l=l, l_N = lN, h=h, h_N=hN, c=c, c_N=c)
    
    # Convert to single-shooting problem
    ocp = OCProblemData(f=f, l=l, l_N=lN, h=h, h_N=hN, c=c, N=N)
    ss_cost, ss_constr = convert_to_single_shooting(ocp)
    cs_u = ss_cost.sx_in(0)

    # Random initialization
    u = rng.uniform(-2, 2, (N * nu,))
    y = rng.uniform(-0.5, 0.5, (nx * N + nx,))
    μ = rng.exponential(1, (nx * N + nx,))
    problem.x_init = rng.uniform(-1, 1, (nx,))

    # Check cost and gradient
    eval = pa.OCPEvaluator(problem)
    V, grad = eval.forward_backward(u, y, μ)

    cs_V = ss_cost(u, problem.x_init)
    assert abs(V - cs_V) < tol

    cs_grad = cs.Function("gr",
        [cs_u],
        [cs.gradient(ss_cost(cs_u, problem.x_init), cs_u)])(u)
    assert la.norm(grad - cs_grad) < tol

    # Check cost and gradient with constraints
    problem.U.lowerbound = rng.uniform(-2.5, 0, nu)
    problem.U.upperbound = rng.uniform(0, +2.5, nu)
    problem.D.lowerbound = rng.uniform(-1, -0.5, nx)
    problem.D.upperbound = rng.uniform(+0.5, +1, nx)
    problem.D_N.lowerbound = problem.D.lowerbound[::-1]
    problem.D_N.upperbound = problem.D.upperbound[::-1]

    aug_eval = pa.OCPEvaluator(problem)
    aug_V, aug_grad = aug_eval.forward_backward(u, y, μ)

    ss_D_lb = np.concatenate((np.tile(problem.D.lowerbound, N), problem.D_N.lowerbound))
    ss_D_ub = np.concatenate((np.tile(problem.D.upperbound, N), problem.D_N.upperbound))
    ss_μ = μ.copy()
    proj_D = lambda g: cs.fmax(ss_D_lb, cs.fmin(g, ss_D_ub))
    ζ = ss_constr(cs_u, problem.x_init) + y / ss_μ
    cs_aug_lagr = cs.Function("al", [cs_u],
        [ss_cost(cs_u, problem.x_init) + 
         0.5 * cs.sum1(ss_μ * (ζ - proj_D(ζ))**2)])

    cs_aug_V = cs_aug_lagr(u)
    assert abs(aug_V - cs_aug_V) < tol

    cs_aug_grad = cs.Function("algr", [cs_u],
        [cs.gradient(cs_aug_lagr(cs_u), cs_u)])(u).full().ravel()
    print("‖∇ψ(C++) - ∇ψ(CasADi) =", la.norm(aug_grad - cs_aug_grad))
    assert la.norm(aug_grad - cs_aug_grad) < tol

    # Simulate
    nxu = nx + nu
    xu = np.nan * np.ones((nxu * N + nx,))
    xu[0:nx] = problem.x_init
    for i in range(N):
        xi = xu[i * nxu:i * nxu + nx]
        ui = u[i * nu:i * nu + nu]
        xu[i * nxu + nx:i * nxu + nxu] = ui
        xu[(i + 1) * nxu:(i + 1) * nxu + nx] = f(xi, ui).full().ravel()

    # Build QP
    n_qp = nxu * N + nx
    Q_qp = np.zeros((n_qp, n_qp))
    for i in range(N):
        Qk = aug_eval.Qk(i, u, y, μ)
        Rk = aug_eval.Rk(i, u, np.arange(nu))
        Sk = aug_eval.Sk(i, u, np.arange(nu))
        Q_qp[i * nxu:i * nxu + nx, i * nxu:i * nxu + nx] = Qk
        Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu + nx:i * nxu + nxu] = Rk
        Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu:i * nxu + nx] = Sk
        Q_qp[i * nxu:i * nxu + nx, i * nxu + nx:i * nxu + nxu] = Sk.T
    Q_qp[N * nxu:N * nxu + nx, N * nxu:N * nxu + nx] = aug_eval.Qk(N, u, y, μ)

    ms_cost_alm = convert_to_multiple_shooting_alm(ocp, y, μ, problem.D, problem.D_N)
    cs_xu = ms_cost_alm.sx_in(0)
    # Compute true Hessian
    cs_hess_qp, _ = cs.hessian(ms_cost_alm(cs_xu), cs_xu)
    cs_hess_qp = cs.evalf(cs.substitute(cs_hess_qp, cs_xu, xu)).full()

    # Gauss-Newton Hessian approximation
    cs_Q_qp = np.zeros((n_qp, n_qp))
    cs_qr_qp = np.nan * np.ones((N * nxu + nx,))
    for i in range(N):
        xi = xu[i * nxu:i * nxu + nx]
        ui = xu[i * nxu + nx:i * nxu + nxu]
        μi = μ[i * nc:i * nc + nc]
        Jhx = cs.substitute(cs.jacobian(h(x_sym, ui), x_sym), x_sym, xi)
        Jhu = cs.substitute(cs.jacobian(h(xi, u_sym), u_sym), u_sym, ui)
        Lhh = cs.substitute(cs.hessian(l(hsym), hsym)[0], hsym, h(xi, ui))
        proj_D = lambda g: cs.fmax(problem.D.lowerbound, cs.fmin(g, problem.D.upperbound))
        c_sym = cs.SX.sym("c", nx)
        ζ = c_sym + y[i * nx:i * nx + nx] / μi
        penalty = 0.5 * cs.sum1(μi * (ζ - proj_D(ζ))**2)
        ci = ocp.c(xi)
        M = cs.substitute(cs.hessian(penalty, c_sym)[0], c_sym, ci)
        Jc = cs.substitute(cs.jacobian(ocp.c(x_sym), x_sym), x_sym, xi)
        Qk = Jhx.T @ Lhh @ Jhx + Jc.T @ M @ Jc
        Rk = Jhu.T @ Lhh @ Jhu
        Sk = Jhu.T @ Lhh @ Jhx
        cs_Q_qp[i * nxu:i * nxu + nx, i * nxu:i * nxu + nx] = cs.evalf(Qk)
        cs_Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu + nx:i * nxu + nxu] = cs.evalf(Rk)
        cs_Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu:i * nxu + nx] = cs.evalf(Sk)
        cs_Q_qp[i * nxu:i * nxu + nx, i * nxu + nx:i * nxu + nxu] = cs.evalf(Sk).T
        grad_l = cs.substitute(cs.gradient(l(hsym), hsym), hsym, h(xi, ui))
        ζ = ci + y[i * nx:i * nx + nx] / μi
        qi = Jhx.T @ grad_l + Jc.T @ (μi * (ζ - proj_D(ζ)))
        ri = Jhu.T @ grad_l
        cs_qr_qp[i * nxu:i * nxu + nx] = cs.evalf(qi).full().ravel()
        cs_qr_qp[i * nxu + nx:i * nxu + nxu] = cs.evalf(ri).full().ravel()
    μN = μ[N * nc:N * nc + nc_N]
    xN = xu[N * nxu:N * nxu + nx]
    JhN = cs.substitute(cs.jacobian(hN(x_sym), x_sym), x_sym, xN)
    LNhh = cs.substitute(cs.hessian(lN(hNsym), hNsym)[0], hNsym, hN(xN))
    proj_D_N = lambda g: cs.fmax(problem.D_N.lowerbound, cs.fmin(g, problem.D_N.upperbound))
    c_sym = cs.SX.sym("c", nx) # TODO: terminal constraints
    ζ = c_sym + y[N * nx:N * nx + nx] / μN
    penalty = 0.5 * cs.sum1(μN * (ζ - proj_D_N(ζ))**2)
    cN = ocp.c(xN)
    MN = cs.substitute(cs.hessian(penalty, c_sym)[0], c_sym, cN)
    JcN = cs.substitute(cs.jacobian(ocp.c(x_sym), x_sym), x_sym, xN)
    QN = JhN.T @ LNhh @ JhN + JcN.T @ MN @ JcN
    cs_Q_qp[N * nxu:N * nxu + nx, N * nxu:N * nxu + nx] = cs.evalf(QN)
    grad_lN = cs.substitute(cs.gradient(lN(hNsym), hNsym), hNsym, hN(xN))
    ζ = cN + y[N * nx:N * nx + nx] / μN
    qN = JhN.T @ grad_lN + JcN.T @ (μN * (ζ - proj_D_N(ζ)))
    cs_qr_qp[N * nxu:N * nxu + nx] = cs.evalf(qN).full().ravel()

    print("‖Q(C++) - Q(CasADi)‖ =", la.norm(Q_qp - cs_Q_qp))
    if "PYTEST_CURRENT_TEST" not in os.environ:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.imshow(np.log10(abs(Q_qp - cs_Q_qp)))
        plt.colorbar()
        plt.title("C++ vs CasADi GN Hessian")
        plt.show()
    assert la.norm(Q_qp - cs_Q_qp) < tol

    if "PYTEST_CURRENT_TEST" not in os.environ:
        plt.figure()
        plt.imshow(np.log10(abs(cs_hess_qp - cs_Q_qp)))
        plt.colorbar()
        plt.title("GN Hessian vs True Hessian")
        plt.show()

    # Build constraint matrix
    triplets = []
    def add_identity(r, c, n, fac, triplets):
        for i in range(n):
            triplets += [(r + i, c + i, fac)]
    def add_block(r, c, M, triplets):
        for col in range(M.shape[1]):
            for row in range(M.shape[0]):
                val = M[row, col]
                if val:
                    triplets += [(r + row, c + col, val)]
    # Dynamics
    for i in range(N + 1):
        add_identity(i * nx, i * nxu, nx, 1., triplets)
        if i > 0:
            j = i - 1
            xj = xu[j * nxu:j * nxu + nx]
            uj = xu[j * nxu + nx:j * nxu + nxu]
            Aj = cs.substitute(cs.jacobian(f(x_sym, uj), x_sym), x_sym, xj)
            Bj = cs.substitute(cs.jacobian(f(xj, u_sym), u_sym), u_sym, uj)
            add_block(i * nx, j * nxu, cs.evalf(-Aj).full(), triplets)
            add_block(i * nx, j * nxu + nx, cs.evalf(-Bj).full(), triplets)
    # Box constraints
    γ = 0.56789
    b = []
    m = 0
    for i in range(N):
        for j in range(nu):
            ui = u[i * nu + j]
            lb_j = problem.U.lowerbound[j]
            ub_j = problem.U.upperbound[j]
            gs = ui - γ * aug_grad[i * nu + j]
            if gs <= lb_j:
                triplets += [(N * nx + nx + m, i * nxu + nx + j, 1.)]
                b += [lb_j - ui]
                m += 1
            elif gs >= ub_j:
                triplets += [(N * nx + nx + m, i * nxu + nx + j, 1.)]
                b += [ub_j - ui]
                m += 1
    triplets = np.array(triplets)
    A_qp = sp.csc_matrix((triplets[:, 2], (triplets[:, 0], triplets[:, 1])), shape=(N * nx + nx + m, N * nxu + nx))

    if "PYTEST_CURRENT_TEST" not in os.environ:
        plt.figure()
        plt.imshow(np.log10(abs(A_qp.toarray())))
        plt.colorbar()
        plt.title("Constraint matrix A")
        plt.show()

    # Solve the QP
    import qpalm
    qp_data = qpalm.Data(Q_qp.shape[0], A_qp.shape[0])
    qp_data.Q, qp_data.A, qp_data.q, qp_data.c = Q_qp, A_qp, cs_qr_qp, 0
    qp_data.bmin = np.concatenate((np.zeros((N * nx + nx,)), np.array(b)))
    qp_data.bmax = np.concatenate((np.zeros((N * nx + nx,)), np.array(b)))

    qp_settings = qpalm.Settings()
    qp_settings.eps_abs = tol * 1e-2
    qp_settings.eps_abs_in = tol * 1e-2
    qp_settings.eps_rel = tol * 1e-2
    qp_settings.eps_rel_in = tol * 1e-2
    qp_settings.eps_prim_inf = tol * 1e-2
    qp_settings.eps_dual_inf = tol * 1e-2
    qp_settings.enable_dual_termination = True
    qp_sol = qpalm.Solver(qp_data, qp_settings)
    qp_sol.solve()
    qp_delxu = qp_sol.solution.x
    assert la.norm(A_qp @ qp_delxu - qp_data.bmin) < tol
    qp_delu = np.reshape(qp_delxu[:-nx], (nxu, N), order='F')[nx:, :]
    assert la.norm(qp_delxu[0:nx]) < tol
    A0 = cs.evalf(cs.substitute(cs.jacobian(f(x_sym, u[:nu]), x_sym), x_sym, xu[:nx])).full()
    B0 = cs.evalf(cs.substitute(cs.jacobian(f(xu[:nx], u_sym), u_sym), u_sym, u[:nu])).full()
    assert la.norm(A0 @ qp_delxu[0:nx] + B0 @ qp_delxu[nx:nx+nu] - qp_delxu[nxu:nxu+nx]) < tol

    # Compare with the LQR solution in C++
    delu = np.reshape(aug_eval.lqr_factor_solve(u, γ, y, μ), (nu, N), order='F')

    # Manual LQR solution (does not handle constraints, unused)
    qN = cs_qr_qp[N * nxu:N * nxu + nx]
    QN = cs_Q_qp[N * nxu:N * nxu + nx, N * nxu:N * nxu + nx]
    s = qN
    P = QN
    K = [None] * N
    e = [None] * N
    A = [None] * N
    B = [None] * N
    for i in range(N - 1, -1, -1):
        xi = xu[i * nxu:i * nxu + nx]
        ui = xu[i * nxu + nx:i * nxu + nxu]
        qi = cs_qr_qp[i * nxu:i * nxu + nx]
        ri = cs_qr_qp[i * nxu + nx:i * nxu + nxu]
        Qi = cs_Q_qp[i * nxu:i * nxu + nx, i * nxu:i * nxu + nx]
        Ri = cs_Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu + nx:i * nxu + nxu]
        Si = cs_Q_qp[i * nxu + nx:i * nxu + nxu, i * nxu:i * nxu + nx]
        Ai = cs.evalf(cs.substitute(cs.jacobian(f(x_sym, ui), x_sym), x_sym, xi)).full()
        Bi = cs.evalf(cs.substitute(cs.jacobian(f(xi, u_sym), u_sym), u_sym, ui)).full()
        R̅ = Ri + Bi.T @ P @ Bi
        S̅ = Si + Bi.T @ P @ Ai
        y = s
        t = Bi.T @ y + ri
        K[i] = -la.solve(R̅, S̅)
        e[i] = -la.solve(R̅, t)
        # print("R̅", R̅, "S̅", S̅, "y", y, "t", t, "Ki", K[i], "ei", e[i], sep='\n')
        s = S̅.T @ e[i] + Ai.T @ y + qi
        P = Qi + Ai.T @ P @ Ai + S̅.T @ K[i]
        A[i] = Ai
        B[i] = Bi
    py_delxu = np.zeros((N * nxu + nx,))
    for i in range(N):
        delxi = py_delxu[i * nxu:i * nxu + nx]
        delui = py_delxu[i * nxu + nx:i * nxu + nxu]
        delui[:] = K[i] @ delxi + e[i]
        py_delxu[(i+1) * nxu:(i+1) * nxu + nx] = A[i] @ delxi + B[i] @ delui

    # Show difference between QP and LQR solutions
    qp_delU = qp_delu.reshape(-1, order='F')
    delU = delu.reshape(-1, order='F')
    print("‖Δu(C++) - Δu(QPALM)‖ =", la.norm(delU - qp_delU))
    if "PYTEST_CURRENT_TEST" not in os.environ:
        plt.figure()
        plt.plot(qp_delU, '.-', label="QPALM")
        plt.plot(delU, '.-', label="C++ LQR")
        plt.title("QPALM solution vs C++ LQR solution")
        plt.legend()
        plt.show()
    assert la.norm(delU - qp_delU) < 1e3 * tol

if __name__ == '__main__':
    test_lqr()
