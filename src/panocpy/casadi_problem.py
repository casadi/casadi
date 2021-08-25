from typing import Tuple
import casadi as cs
from tempfile import TemporaryDirectory
import os
from panocpy import load_casadi_problem_with_param
from panocpy import Problem


def generate_casadi_problem(
    name: str, f: cs.Function, g: cs.Function, second_order: bool = False
) -> Tuple[cs.CodeGenerator, int, int, int]:

    assert f.n_in() in [1, 2]
    assert f.n_in() == g.n_in()
    assert f.size1_in(0) == g.size1_in(0)
    if f.n_in() == 2:
        assert f.size1_in(1) == g.size1_in(1)
    assert f.n_out() == 1
    assert g.n_out() == 1
    n = f.size1_in(0)
    m = g.size1_out(0)
    p = f.size1_in(1) if f.n_in() == 2 else 0
    print(f"{n=}, {m=}, {p=}")
    xp = (f.sx_in(0), f.sx_in(1)) if f.n_in() == 2 else (f.sx_in(0),)
    xp_names = (f.name_in(0), f.name_in(1)) if f.n_in() == 2 else (f.name_in(0),)
    x = xp[0]
    y = cs.SX.sym("y", m)
    v = cs.SX.sym("v", n)

    L = f(*xp) + cs.dot(y, g(*xp)) if m > 0 else f(*xp)

    cgname = f"{name}.c"
    cg = cs.CodeGenerator(cgname)
    cg.add(
        cs.Function(
            "f",
            [*xp],
            [f(*xp)],
            [*xp_names],
            ["f"],
        )
    )
    cg.add(
        cs.Function(
            "grad_f",
            [*xp],
            [cs.gradient(f(*xp), x)],
            [*xp_names],
            ["grad_f"],
        )
    )
    cg.add(
        cs.Function(
            "g",
            [*xp],
            [g(*xp)],
            [*xp_names],
            ["g"],
        )
    )
    cg.add(
        cs.Function(
            "grad_g",
            [*xp, y],
            [cs.jtimes(g(*xp), x, y, True)],
            [*xp_names, "y"],
            ["grad_g"],
        )
    )
    if second_order:
        cg.add(
            cs.Function(
                "hess_L",
                [*xp, y],
                [cs.hessian(L, x)[0]],
                [*xp_names, "y"],
                ["hess_L"],
            )
        )
        cg.add(
            cs.Function(
                "hess_L_prod",
                [*xp, y, v],
                [cs.gradient(cs.jtimes(L, x, v, False), x)],
                [*xp_names, "y", "v"],
                ["hess_L_prod"],
            )
        )
    return cg, n, m, p


def compile_and_load_problem(
    cgen: cs.CodeGenerator,
    n: int,
    m: int,
    name: str = "PANOC_ALM_problem",
) -> Problem:
    """Compile the C-code using the given code-generator and load it as a
    panocpy Problem.

    Args:
        cgen (cs.CodeGenerator): Code generator to generate C-code for the costs and the constraints with.
        n (int): Dimensions of the decision variables (primal dimension)
        m (int): Number of nonlinear constraints (dual dimension)
        name (str, optional): String description of the problem. Defaults to "PANOC_ALM_problem".
    """

    with TemporaryDirectory(prefix="") as tmpdir:
        cfile = cgen.generate(tmpdir)
        sofile = os.path.join(tmpdir, f"{name}.so")
        os.system(f"cc -fPIC -shared -O3 -march=native {cfile} -o {sofile}")
        prob = load_casadi_problem_with_param(sofile, n, m)
    return prob
