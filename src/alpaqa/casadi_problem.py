from typing import Tuple, Union
import casadi as cs
import alpaqa as pa
import os
from os.path import join, basename
import shelve
import uuid
import pickle
import base64
import glob
import subprocess
import tempfile
import platform


def generate_casadi_problem(
    f: cs.Function,
    g: cs.Function,
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
    assert f.n_in() == g.n_in()
    assert f.size1_in(0) == g.size1_in(0)
    if f.n_in() == 2:
        assert f.size1_in(1) == g.size1_in(1)
    assert f.n_out() == 1
    assert g.n_out() == 1
    n = f.size1_in(0)
    m = g.size1_out(0)
    p = f.size1_in(1) if f.n_in() == 2 else 0
    xp = (f.sx_in(0), f.sx_in(1)) if f.n_in() == 2 else (f.sx_in(0), )
    xp_names = (f.name_in(0),
                f.name_in(1)) if f.n_in() == 2 else (f.name_in(0), )
    x = xp[0]
    y = cs.SX.sym("y", m)
    v = cs.SX.sym("v", n)

    L = f(*xp) + cs.dot(y, g(*xp)) if m > 0 else f(*xp)

    cgname = f"{name}.c"
    cg = cs.CodeGenerator(cgname)
    cg.add(cs.Function(
        "f",
        [*xp],
        [f(*xp)],
        [*xp_names],
        ["f"],
    ))
    cg.add(
        cs.Function(
            "grad_f",
            [*xp],
            [cs.gradient(f(*xp), x)],
            [*xp_names],
            ["grad_f"],
        ))
    cg.add(cs.Function(
        "g",
        [*xp],
        [g(*xp)],
        [*xp_names],
        ["g"],
    ))
    cg.add(
        cs.Function(
            "grad_g",
            [*xp, y],
            [cs.jtimes(g(*xp), x, y, True)],
            [*xp_names, "y"],
            ["grad_g"],
        ))
    if second_order:
        cg.add(
            cs.Function(
                "hess_L",
                [*xp, y],
                [cs.hessian(L, x)[0]],
                [*xp_names, "y"],
                ["hess_L"],
            ))
        cg.add(
            cs.Function(
                "hess_L_prod",
                [*xp, y, v],
                [cs.gradient(cs.jtimes(L, x, v, False), x)],
                [*xp_names, "y", "v"],
                ["hess_L_prod"],
            ))
    return cg, n, m, p


def _load_casadi_problem(sofile, n, m, p):
    if p > 0:
        prob = pa.load_casadi_problem_with_param(sofile, n, m, p)
    else:
        prob = pa.load_casadi_problem(sofile, n, m)
    return prob


def generate_and_compile_casadi_problem(
    f: cs.Function,
    g: cs.Function,
    second_order: bool = False,
    name: str = "alpaqa_problem",
) -> Union[pa.Problem, pa.ProblemWithParam]:
    """Compile the objective and constraint functions into a alpaqa Problem.

    :param f:            Objective function.
    :param g:            Constraint function.
    :param second_order: Whether to generate functions for evaluating Hessians.
    :param name: Optional string description of the problem (used for filename).

    :return:   * Problem specification that can be passed to the solvers.
    """

    cachedir = None
    if not cachedir:
        homecachedir = os.path.expanduser("~/.cache")
        if os.path.isdir(homecachedir):
            cachedir = join(homecachedir, 'alpaqa', 'cache')
        else:
            cachedir = join(tempfile.gettempdir(), 'alpaqa', 'cache')
    cachefile = join(cachedir, 'problems')

    key = base64.b64encode(pickle.dumps(
        (f, g, second_order, name))).decode('ascii')

    os.makedirs(cachedir, exist_ok=True)
    with shelve.open(cachefile) as cache:
        if key in cache:
            uid, soname, dim = cache[key]
            probdir = join(cachedir, str(uid))
            sofile = join(probdir, soname)
            try:
                return _load_casadi_problem(sofile, *dim)
            except:
                del cache[key]
                # if os.path.exists(probdir) and os.path.isdir(probdir):
                #     shutil.rmtree(probdir)
                raise
        uid = uuid.uuid1()
        projdir = join(cachedir, "build")
        builddir = join(projdir, "build")
        os.makedirs(builddir, exist_ok=True)
        probdir = join(cachedir, str(uid))
        cgen, n, m, p = generate_casadi_problem(f, g, second_order, name)
        cfile = cgen.generate(join(projdir, ""))
        with open(join(projdir, 'CMakeLists.txt'), 'w') as f:
            f.write(f"""
                cmake_minimum_required(VERSION 3.17)
                project(CasADi-{name} LANGUAGES C)
                set(CMAKE_SHARED_LIBRARY_PREFIX "")
                add_library({name} SHARED {basename(cfile)})
                install(FILES $<TARGET_FILE:{name}>
                        DESTINATION lib)
                install(FILES {basename(cfile)}
                        DESTINATION src)
            """)
        build_type = 'Release'
        configure_cmd = ['cmake', '-B', builddir, '-S', projdir]
        if platform.system() != 'Windows':
            configure_cmd += ['-G', 'Ninja Multi-Config']
        build_cmd = ['cmake', '--build', builddir, '--config', build_type]
        install_cmd = [
            'cmake', '--install', builddir, '--config', build_type, '--prefix',
            probdir
        ]
        subprocess.run(configure_cmd, check=True)
        subprocess.run(build_cmd, check=True)
        subprocess.run(install_cmd, check=True)
        sofile = glob.glob(join(probdir, "lib", name + ".*"))
        if len(sofile) == 0:
            raise RuntimeError(
                f"Unable to find compiled CasADi problem '{name}'")
        elif len(sofile) > 1:
            raise RuntimeWarning(
                f"Multiple compiled CasADi problem files were found for '{name}'"
            )
        sofile = sofile[0]
        soname = os.path.relpath(sofile, probdir)
        cache[key] = uid, soname, (n, m, p)

        return _load_casadi_problem(sofile, n, m, p)