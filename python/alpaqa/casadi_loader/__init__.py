import casadi as cs
import os
from os.path import join, basename
import alpaqa as pa
import shelve
import uuid
import pickle
import base64
import glob
import subprocess
import tempfile
import platform
from ..casadi_generator import generate_casadi_problem


def _load_casadi_problem(sofile, n, m, p):
    print("-- Loading:", sofile)
    prob = pa.load_casadi_problem(sofile, n=n, m=m, p=p)
    return prob


def generate_and_compile_casadi_problem(
    f: cs.Function,
    g: cs.Function,
    second_order: bool = False,
    name: str = "alpaqa_problem",
) -> pa.CasADiProblem:
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
        if platform.system() == 'Windows':
            configure_cmd += ['-A', 'x64' if sys.maxsize > 2**32 else 'Win32']
        else:
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