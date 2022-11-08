import casadi as cs
from typing import Tuple, Optional


def generate_casadi_problem(
    f: cs.Function,
    g: Optional[cs.Function],
    second_order: bool = False,
    name: str = "alpaqa_problem",
) -> Tuple[cs.CodeGenerator, int, int, int]:
    """Convert the objective and constraint functions into a CasADi code
    generator.

    :param f:            Objective function.
    :param g:            Constraint function.
    :param second_order: Whether to generate functions for evaluating Hessians.
    :param name: Optional string description of the problem (used for filename).

    :return:   * Code generator that generates the functions and derivatives
                 used by the solvers.
               * Dimensions of the decision variables (primal dimension).
               * Number of nonlinear constraints (dual dimension).
               * Number of parameters.
    """

    assert f.n_in() in [1, 2]
    assert f.n_out() == 1
    n = f.size1_in(0)
    if g is not None:
        assert f.n_in() == g.n_in()
        assert f.size1_in(0) == g.size1_in(0)
        if f.n_in() == 2:
            assert f.size1_in(1) == g.size1_in(1)
        assert g.n_out() <= 1
        m = g.size1_out(0) if g.n_out() == 1 else 0
    else:
        m = 0
    p = f.size1_in(1) if f.n_in() == 2 else 0
    xp = (f.sx_in(0), f.sx_in(1)) if f.n_in() == 2 else (f.sx_in(0), )
    xp_def = (f.sx_in(0), f.sx_in(1)) if f.n_in() == 2 \
        else (f.sx_in(0), cs.SX.sym("p", 0))
    xp_names = (f.name_in(0), f.name_in(1)) if f.n_in() == 2 \
          else (f.name_in(0), "p")
    x = xp[0]
    y = cs.SX.sym("y", m)
    v = cs.SX.sym("v", n)
    Σ = cs.SX.sym("Σ", m)
    zl = cs.SX.sym("zl", m)
    zu = cs.SX.sym("zu", m)

    if m > 0:
        L = f(*xp) + cs.dot(y, g(*xp))
        ζ = g(*xp) + (y / Σ)
        ẑ = cs.fmax(zl, cs.fmin(ζ, zu))
        d = ζ - ẑ
        ŷ = Σ * d
        ψ = f(*xp) + 0.5 * cs.dot(ŷ, d)
    else:
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

    if second_order:
        cg.add(
            cs.Function(
                "hess_L",
                [*xp_def, y],
                [cs.hessian(L, x)[0]],
                [*xp_names, "y"],
                ["hess_L"],
            ))
        cg.add(
            cs.Function(
                "hess_L_prod",
                [*xp_def, y, v],
                [cs.gradient(cs.jtimes(L, x, v, False), x)],
                [*xp_names, "y", "v"],
                ["hess_L_prod"],
            ))
    return cg, n, m, p

def generate_casadi_control_problem(
    f: cs.Function,
    l: cs.Function,
    l_N: cs.Function,
    name: str = "alpaqa_control_problem",
) -> Tuple[cs.CodeGenerator, int, int, int]:
    """Convert the dynamics and cost functions into a CasADi code generator.

    :param f:            Dynamics.
    :param l:            Stage cost.
    :param l_N:          Terminal cost.
    :param name: Optional string description of the problem (used for filename).

    :return:   * Code generator that generates the functions and derivatives
                 used by the solvers.
               * Dimensions of the decision variables (primal dimension).
               * Number of nonlinear constraints (dual dimension).
               * Number of parameters.
    """

    assert f.n_in() in [2, 3]
    assert f.n_out() == 1
    nx = f.size1_in(0)
    nu = f.size1_in(1)
    p = f.size1_in(2) if f.n_in() == 3 else 0
    assert f.size1_out(0) == nx
    xup = (f.sx_in(0), f.sx_in(1), f.sx_in(2)) if f.n_in() == 3 \
        else (f.sx_in(0), f.sx_in(1))
    xup_def = (f.sx_in(0), f.sx_in(1), f.sx_in(2)) if f.n_in() == 3 \
        else (f.sx_in(0), f.sx_in(1), cs.SX.sym("p", 0))
    xup_names = (f.name_in(0), f.name_in(1), f.name_in(2)) if f.n_in() == 3 \
          else (f.name_in(0), f.name_in(1), "p")
    x = xup[0]
    u = xup[1]

    nh = nx + nu
    assert l.n_in() == f.n_in() - 1
    l_p = l.size1_in(1) if l.n_in() == 2 else 0
    assert p == l_p
    assert l.size1_in(0) == nh
    assert l.size1_out(0) == 1

    hp = (l.sx_in(0), l.sx_in(1)) if l.n_in() == 2 \
        else (l.sx_in(0), )
    hp_def = (l.sx_in(0), l.sx_in(1)) if l.n_in() == 2 \
        else (l.sx_in(0), cs.SX.sym("p", 0))
    hp_names = (l.name_in(0), l.name_in(1)) if l.n_in() == 2 \
          else (l.name_in(0), "p")
    h = hp[0]

    assert l_N.n_in() == f.n_in() - 1
    l_N_p = l_N.size1_in(1) if l_N.n_in() == 2 else 0
    assert p == l_N_p
    assert l_N.size1_in(0) == nx
    assert l_N.size1_out(0) == 1

    xp = (l_N.sx_in(0), l_N.sx_in(1)) if l_N.n_in() == 2 \
        else (l_N.sx_in(0), )
    xp_def = (l_N.sx_in(0), l_N.sx_in(1)) if l_N.n_in() == 2 \
        else (l_N.sx_in(0), cs.SX.sym("p", 0))
    xp_names = (l_N.name_in(0), l_N.name_in(1)) if l_N.n_in() == 2 \
          else (l_N.name_in(0), "p")
    lNx = xp[0]

    v = cs.SX.sym("v", nx)

    cgname = f"{name}.c"
    cg = cs.CodeGenerator(cgname)
    cg.add(cs.Function(
        "f",
        [*xup_def],
        [f(*xup)],
        [*xup_names],
        ["f"],
    ))
    cg.add(cs.Function(
        "jac_f",
        [*xup_def],
        [cs.densify(cs.jacobian(f(*xup), cs.vertcat(x, u)))],
        [*xup_names],
        ["jac_f"],
    ))
    cg.add(cs.Function(
        "grad_f_prod",
        [*xup_def, v],
        [cs.jtimes(f(*xup), cs.vertcat(x, u), v, True)],
        [*xup_names, "v"],
        ["grad_f_prod"],
    ))
    cg.add(cs.Function(
        "l",
        [*hp_def],
        [l(*hp)],
        [*hp_names],
        ["l"],
    ))
    cg.add(cs.Function(
        "l_N",
        [*xp_def],
        [l_N(*xp)],
        [*xp_names],
        ["l_N"],
    ))
    cg.add(cs.Function(
        "grad_l",
        [*hp_def],
        [cs.gradient(l(*hp), h)],
        [*hp_names],
        ["grad_l"],
    ))
    cg.add(cs.Function(
        "grad_l_N",
        [*xp_def],
        [cs.gradient(l_N(*xp), lNx)],
        [*xp_names],
        ["grad_l_N"],
    ))
    cg.add(cs.Function(
        "hess_l",
        [*hp_def],
        [cs.densify(cs.hessian(l(*hp), h)[0])],
        [*hp_names],
        ["hess_l"],
    ))
    cg.add(cs.Function(
        "hess_l_N",
        [*xp_def],
        [cs.densify(cs.hessian(l_N(*xp), lNx)[0])],
        [*xp_names],
        ["hess_l_N"],
    ))
    return cg, nx, nu, p