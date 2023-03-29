import casadi as cs
import numpy as np
from os.path import splitext
from typing import Tuple, Optional, Literal, get_args, Callable

SECOND_ORDER_SPEC = Literal['no', 'full', 'prod', 'L', 'psi_prod', 'psi']

def generate_casadi_problem(
    f: cs.Function,
    g: Optional[cs.Function],
    second_order: SECOND_ORDER_SPEC = 'no',
    name: str = "alpaqa_problem",
    sym: Callable = cs.SX.sym,
) -> cs.CodeGenerator:
    """Convert the objective and constraint functions into a CasADi code
    generator.

    :param f:            Objective function.
    :param g:            Constraint function.
    :param second_order: Whether to generate functions for evaluating Hessians.
    :param name: Optional string description of the problem (used for filename).

    :return: Code generator that generates the functions and derivatives used
             by the solvers.
    """

    assert second_order in get_args(SECOND_ORDER_SPEC)

    assert f.n_in() in [1, 2]
    assert f.n_out() == 1
    assert f.size1_out(0) == 1
    assert f.size2_out(0) == 1
    n = f.size1_in(0)
    assert f.size2_in(0) == 1
    with_param = f.n_in() == 2
    if with_param:
        assert f.size2_in(1) == 1
    if g is not None:
        assert f.n_in() == g.n_in()
        assert f.size1_in(0) == g.size1_in(0)
        if with_param:
            assert f.size1_in(1) == g.size1_in(1)
            assert g.size2_in(1) == 1
        assert g.n_out() <= 1
        m = g.size1_out(0) if g.n_out() == 1 else 0
        if g.n_out() == 1:
            assert g.size2_out(0) == 1
    else:
        m = 0
    x = sym("x", n)
    p = sym("p", f.size1_in(1) if with_param else 0)
    xp = (x, p) if with_param else (x, )
    xp_def = (x, p)
    xp_names = (f.name_in(0), f.name_in(1)) if with_param \
          else (f.name_in(0), "p")
    x = xp[0]
    y = sym("y", m)
    s = sym("s")
    v = sym("v", n)
    Σ = sym("Σ", m)
    zl = sym("zl", m)
    zu = sym("zu", m)

    if m > 0:
        sL = s * f(*xp) + cs.dot(y, g(*xp))
        L = f(*xp) + cs.dot(y, g(*xp))
        ζ = g(*xp) + (y / Σ)
        ẑ = cs.fmax(zl, cs.fmin(ζ, zu))
        d = ζ - ẑ
        ŷ = Σ * d
        sψ = s * f(*xp) + 0.5 * cs.dot(ŷ, d)
        ψ = f(*xp) + 0.5 * cs.dot(ŷ, d)
    else:
        sL = s * f(*xp)
        sψ = s * f(*xp)
        L = f(*xp)
        ψ = f(*xp)

    cgname = f"{name}.c"
    cg = cs.CodeGenerator(cgname)
    cg.add(cs.Function(
        "f",
        [*xp_def],
        [f(*xp)],
        [*xp_names],
        ["f"],
    ))
    cg.add(
        cs.Function(
            "f_grad_f",
            [*xp_def],
            [f(*xp), cs.gradient(f(*xp), x)],
            [*xp_names],
            ["f", "grad_f"],
        ))
    cg.add(
        cs.Function(
            "g",
            [*xp_def],
            [g(*xp)] if m > 0 else [],
            [*xp_names],
            ["g"] if m > 0 else [],
        ))
    cg.add(
        cs.Function(
            "psi_grad_psi",
            [*xp_def, y, Σ, zl, zu],
            [ψ, cs.gradient(ψ, x)],
            [*xp_names, "y", "Σ", "zl", "zu"],
            ["ψ", "grad_ψ"],
        ))
    if False:
        cg.add(
            cs.Function(
                "grad_psi",
                [*xp_def, y, Σ, zl, zu],
                [cs.gradient(ψ, x)],
                [*xp_names, "y", "Σ", "zl", "zu"],
                ["grad_ψ"],
            ))
    if m > 0:
        cg.add(
            cs.Function(
                "grad_L",
                [*xp_def, y],
                [cs.gradient(L, x)],
                [*xp_names, "y"],
                ["grad_L"],
            ))
        cg.add(
            cs.Function(
                "psi",
                [*xp_def, y, Σ, zl, zu],
                [ψ, ŷ],
                [*xp_names, "y", "Σ", "zl", "zu"],
                ["ψ", "ŷ"],
            ))

    if second_order in ['full', 'L']:
        cg.add(
            cs.Function(
                "jacobian_g",
                [*xp_def],
                [cs.jacobian(g(*xp), x)],
                [*xp_names],
                ["jac_g"],
            ))
        HL = cs.hessian(sL, x)[0]
        if not HL.is_dense():
            HL = cs.triu(HL)
        cg.add(
            cs.Function(
                "hess_L",
                [*xp_def, y, s],
                [HL],
                [*xp_names, "y", "s"],
                ["hess_L"],
            ))
    if second_order in ['full', 'prod', 'L']:
        cg.add(
            cs.Function(
                "hess_L_prod",
                [*xp_def, y, s, v],
                [cs.gradient(cs.jtimes(sL, x, v, False), x)],
                [*xp_names, "y", "s", "v"],
                ["hess_L_prod"],
            ))
    if second_order in ['full', 'psi']:
        Hψ = cs.hessian(sψ, x)[0]
        if not Hψ.is_dense():
            Hψ = cs.triu(Hψ)
        cg.add(
            cs.Function(
                "hess_psi",
                [*xp_def, y, Σ, s, zl, zu],
                [Hψ],
                [*xp_names, "y", "Σ", "s", "zl", "zu"],
                ["hess_psi"],
            )
        )
    if second_order in ['full', 'prod', 'psi_prod', 'psi_prod']:
        cg.add(
            cs.Function(
                "hess_psi_prod",
                [*xp_def, y, Σ, s, zl, zu, v],
                [cs.gradient(cs.jtimes(sψ, x, v, False), x)],
                [*xp_names, "y", "Σ", "s", "zl", "zu", "v"],
                ["hess_psi_prod"],
            ))
    return cg

def _add_parameter(f: cs.Function, expected_inputs: int) -> Tuple[cs.Function, cs.SX, str]:
    if f.n_in() == expected_inputs + 1:
        # Okay, we already have a parameter argument
        return f, f.sx_in(expected_inputs), f.name_in(expected_inputs)
    elif f.n_in() == expected_inputs:
        # We don't have a parameter argument
        param = cs.SX.sym("p", 0)
        fx = f(*(f.sx_in(i) for i in range(expected_inputs)))
        return cs.Function(
            f.name(),
            [f.sx_in(i) for i in range(expected_inputs)] + [param],
            [fx] if f.n_out() == 1 else fx,
            [f.name_in(i) for i in range(expected_inputs)] + ["p"],
            [f.name_out(i) for i in range(f.n_out())],
        ), param, "p"
    else:
        raise RuntimeError(f"Incorrect number of inputs for {f.name()} "
                           f"(expected {expected_inputs} inputs with optional "
                           f"additional parameter)")

def generate_casadi_control_problem(
    f: cs.Function,
    l: cs.Function,
    l_N: cs.Function,
    h: cs.Function = None,
    h_N: cs.Function = None,
    c: cs.Function = None,
    c_N: cs.Function = None,
    name: str = "alpaqa_control_problem",
) -> cs.CodeGenerator:
    """Convert the dynamics and cost functions into a CasADi code generator.

    :param f:            Dynamics.
    :param name: Optional string description of the problem (used for filename).

    :return:   * Code generator that generates the functions and derivatives
                 used by the solvers.
               * Number of states.
               * Number of inputs.
               * Number of parameters.
    """

    cgname = f"{name}.c"
    cg = cs.CodeGenerator(cgname)

    assert f.n_in() in [2, 3]
    assert f.n_out() == 1
    f, p_var, p_name = _add_parameter(f, 2) # x, u
    assert f.size2_in(0) == 1
    nx = f.size1_in(0)
    assert f.size2_in(1) == 1
    nu = f.size1_in(1)
    assert f.size2_in(2) == 1
    p = p_var.size1()
    assert f.size1_out(0) == nx
    assert f.size2_out(0) == 1
    x_var = f.sx_in(0)
    u_var = f.sx_in(1)
    xu_var = cs.vertcat(x_var, u_var)

    v_var = cs.SX.sym("v", nx)

    # dynamics and their derivatives
    cg.add(cs.Function(
        "f",
        [x_var, u_var, p_var],
        [f(x_var, u_var, p_var)],
        [f.name_in(i) for i in range(3)],
        [f.name_out(0)],
    ))
    cg.add(cs.Function(
        "jacobian_f",
        [x_var, u_var, p_var],
        [cs.densify(cs.jacobian(f(x_var, u_var, p_var), xu_var))],
        [f.name_in(i) for i in range(3)],
        ["jac_" + f.name_out(0)],
    ))
    cg.add(cs.Function(
        "grad_f_prod",
        [x_var, u_var, p_var, v_var],
        [cs.jtimes(f(x_var, u_var, p_var), xu_var, v_var, True)],
        [f.name_in(i) for i in range(3)] + ["v"],
        ["grad_" + f.name_out(0) + "_prod"],
    ))

    # output mapping
    if h is None:
        h = cs.Function("h", [x_var, u_var, p_var], [xu_var])
    else:
        assert h.n_in() in [2, 3]
        h, _, _ = _add_parameter(h, 2) # x, u
    assert h.size1_in(0) == nx
    assert h.size2_in(0) == 1
    assert h.size1_in(1) == nu
    assert h.size2_in(1) == 1
    assert h.size1_in(2) == p
    assert h.size2_in(2) == 1
    nh = h.size1_out(0)
    assert h.size2_out(0) == 1

    cg.add(cs.Function(
        "h",
        [x_var, u_var, p_var],
        [h(x_var, u_var, p_var)],
        [h.name_in(i) for i in range(3)],
        [h.name_out(0)],
    ))

    # terminal output mapping
    if h_N is None:
        h_N = cs.Function("h_N", [x_var, p_var], [x_var])
    else:
        assert h_N.n_in() in [1, 2]
        h_N, _, _ = _add_parameter(h_N, 1) # x
    assert h_N.size1_in(0) == nx
    assert h_N.size2_in(0) == 1
    assert h_N.size1_in(1) == p
    assert h_N.size2_in(1) == 1
    nh_N = h_N.size1_out(0)
    assert h_N.size2_out(0) == 1

    cg.add(cs.Function(
        "h_N",
        [x_var, p_var],
        [h_N(x_var, p_var)],
        [h_N.name_in(i) for i in range(2)],
        [h_N.name_out(0)],
    ))

    # cost
    assert l.n_in() in [1, 2] # h
    l, _, _ = _add_parameter(l, 1)
    assert l.size1_in(0) == nh
    assert l.size2_in(0) == 1
    assert l.size1_in(1) == p
    assert l.size2_in(1) == 1
    assert l.n_out() == 1
    assert l.size1_out(0) == 1
    assert l.size2_out(0) == 1

    h_var = cs.SX.sym("h", *l.sx_in(0).shape)

    cg.add(cs.Function(
        "l",
        [h_var, p_var],
        [l(h_var, p_var)],
        [l.name_in(i) for i in range(2)],
        [l.name_out(0)],
    ))

    cg.add(cs.Function(
        "qr",
        [xu_var, h_var, p_var],
        [cs.jtimes(h(x_var, u_var, p_var), xu_var, cs.gradient(l(h_var, p_var), h_var), True)],
        ["xu", "h", "p"], # TODO
        ["qr"], # TODO
    ))

    # (JhᵀΛ)ᵀ = ΛJh
    # JhᵀΛ = cs.jtimes(h(x_var, u_var, p_var), x_var, cs.hessian(l(h_var, p_var), h_var)[0], True)
    # JhᵀΛJh = cs.jtimes(h(x_var, u_var, p_var), x_var, cs.transpose(JhTΛ), True)
    Jhx = cs.jacobian(h(x_var, u_var, p_var), x_var)
    Λ = cs.hessian(l(h_var, p_var), h_var)[0]
    Q = Jhx.T @ Λ @ Jhx
    cg.add(cs.Function(
        "Q",
        [xu_var, h_var, p_var],
        [Q],
        ["xu", "h", "p"], # TODO
        ["Q"], # TODO
    ))

    # (JhᵀΛ)ᵀ = ΛJh
    # JhᵀΛ = cs.jtimes(h(x_var, u_var, p_var), u_var, cs.hessian(l(h_var, p_var), h_var)[0], True)
    # JhᵀΛJh = cs.jtimes(h(x_var, u_var, p_var), u_var, cs.transpose(JhTΛ), True)
    Jhu = cs.jacobian(h(x_var, u_var, p_var), u_var)
    R = Jhu.T @ Λ @ Jhu
    cg.add(cs.Function(
        "R",
        [xu_var, h_var, p_var],
        [R],
        ["xu", "h", "p"], # TODO
        ["R"], # TODO
    ))

    # (JhᵀΛ)ᵀ = ΛJh
    # JhᵀΛ = cs.jtimes(h(x_var, u_var, p_var), x_var, cs.hessian(l(h_var, p_var), h_var)[0], True)
    # JhᵀΛJh = cs.jtimes(h(x_var, u_var, p_var), u_var, cs.transpose(JhTΛ), True)
    S = Jhu.T @ Λ @ Jhx
    cg.add(cs.Function(
        "S",
        [xu_var, h_var, p_var],
        [S],
        ["xu", "h", "p"], # TODO
        ["S"], # TODO
    ))

    # terminal cost
    assert l_N.n_in() in [1, 2] # h
    l_N, _, _ = _add_parameter(l_N, 1)
    assert l_N.size1_in(0) == nh_N
    assert l_N.size2_in(0) == 1
    assert l_N.size1_in(1) == p
    assert l_N.size2_in(1) == 1
    assert l_N.n_out() == 1
    assert l_N.size1_out(0) == 1
    assert l_N.size2_out(0) == 1

    hN_var = cs.SX.sym("hN", *l_N.sx_in(0).shape)

    cg.add(cs.Function(
        "l_N",
        [hN_var, p_var],
        [l_N(hN_var, p_var)],
        [l_N.name_in(i) for i in range(2)],
        [l_N.name_out(0)],
    ))

    cg.add(cs.Function(
        "q_N",
        [x_var, hN_var, p_var],
        [cs.jtimes(h_N(x_var, p_var), x_var, cs.gradient(l_N(hN_var, p_var), hN_var), True)],
        ["x", "h", "p"], # TODO
        ["q_N"], # TODO
    ))

    # (JhᵀΛ)ᵀ = ΛJh
    # JhᵀΛ = cs.jtimes(h_N(x_var, p_var), x_var, cs.hessian(l_N(hN_var, p_var), hN_var)[0], True)
    # JhᵀΛJh = cs.jtimes(h_N(x_var, p_var), x_var, cs.transpose(JhTΛ), True)
    JhN = cs.jacobian(h_N(x_var, p_var), x_var)
    ΛN = cs.hessian(l_N(hN_var, p_var), hN_var)[0]
    Q_N = JhN.T @ ΛN @ JhN
    cg.add(cs.Function(
        "Q_N",
        [x_var, hN_var, p_var],
        [Q_N],
        ["x", "h", "p"], # TODO
        ["Q_N"], # TODO
    ))

    # constraints
    if c is None:
        c = cs.Function("c", [x_var, p_var], [cs.vertcat()])
    else:
        assert c.n_in() in [1, 2]
        c, _, _ = _add_parameter(c, 1)
    assert c.size1_in(0) == nx
    assert c.size2_in(0) == 1
    assert c.size1_in(1) == p
    assert c.size2_in(1) == 1
    assert c.n_out() == 1
    nc = c.size1_out(0)
    assert c.size2_out(0) in [0, 1]

    w_var = cs.SX.sym("w", nc)

    cg.add(cs.Function(
        "c",
        [x_var, p_var],
        [c(x_var, p_var)],
        [c.name_in(i) for i in range(2)],
        [c.name_out(0)],
    ))

    cg.add(cs.Function(
        "grad_c_prod",
        [x_var, p_var, w_var],
        [cs.jtimes(c(x_var, p_var), x_var, w_var, True) if nc > 0 else cs.DM.zeros(nx)],
        [c.name_in(i) for i in range(2)] + ["w"],
        ["grad_" + c.name_out(0) + "_prod"],
    ))

    m_var = cs.SX.sym("m", nc)
    # (JhᵀM)ᵀ = MJh
    # JhᵀM = cs.jtimes(c(x_var, p_var), x_var, cs.diag(m_var), True)
    # JhᵀMJh = cs.jtimes(c(x_var, p_var), x_var, cs.transpose(JhᵀM), True)
    Jc = cs.jacobian(c(x_var, p_var), x_var)
    JhᵀMJh = Jc.T @ cs.diag(m_var) @ Jc
    cg.add(cs.Function(
        "gn_hess_c",
        [x_var, p_var, m_var],
        [JhᵀMJh],
        [c.name_in(i) for i in range(2)] + ["m"],
        ["gn_hess_" + c.name_out(0)],
    ))

    # constraints
    if c_N is None:
        c_N = cs.Function("c_N", [x_var, p_var], [cs.vertcat()])
    else:
        assert c_N.n_in() in [1, 2]
        c_N, _, _ = _add_parameter(c_N, 1)
    assert c_N.size1_in(0) == nx
    assert c_N.size2_in(0) == 1
    assert c_N.size1_in(1) == p
    assert c_N.size2_in(1) == 1
    assert c_N.n_out() == 1
    nc_N = c_N.size1_out(0)
    assert c_N.size2_out(0) in [0, 1]

    wN_var = cs.SX.sym("wN", nc_N)

    cg.add(cs.Function(
        "c_N",
        [x_var, p_var],
        [c_N(x_var, p_var)],
        [c_N.name_in(i) for i in range(2)],
        [c_N.name_out(0)],
    ))

    cg.add(cs.Function(
        "grad_c_prod_N",
        [x_var, p_var, wN_var],
        [cs.jtimes(c_N(x_var, p_var), x_var, wN_var, True) if nc_N > 0 else cs.DM.zeros(nx)],
        [c_N.name_in(i) for i in range(2)] + ["w"],
        ["grad_" + c_N.name_out(0) + "_prod"],
    ))

    mN_var = cs.SX.sym("mN", nc_N)
    JcN = cs.jacobian(c_N(x_var, p_var), x_var)
    JhᵀMJhN = JcN.T @ cs.diag(mN_var) @ JcN
    cg.add(cs.Function(
        "gn_hess_c_N",
        [x_var, p_var, mN_var],
        [JhᵀMJhN],
        [c_N.name_in(i) for i in range(2)] + ["m"],
        ["gn_hess_" + c_N.name_out(0)],
    ))

    return cg

def write_casadi_problem_data(sofile, C, D, param):
    if C is None and D is None and param is None:
        return
    C = ([], []) if C is None else C
    D = ([], []) if D is None else D
    param = [] if param is None else param
    with open(f"{splitext(sofile)[0]}.csv", "w") as f:
        opt = dict(delimiter=",", newline="\n")
        ravelrow = lambda x: np.reshape(x, (1, -1), order="A")
        writerow = lambda x: np.savetxt(f, ravelrow(x), **opt)
        writerow(getattr(C, "lowerbound", None) or C[0])
        writerow(getattr(C, "upperbound", None) or C[1])
        writerow(getattr(D, "lowerbound", None) or D[0])
        writerow(getattr(D, "upperbound", None) or D[1])
        writerow(param)
