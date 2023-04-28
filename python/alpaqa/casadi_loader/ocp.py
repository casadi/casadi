import casadi as cs
import os
from os.path import join, basename
import shelve
import uuid
import pickle
import base64
import glob
import subprocess
import platform
import sys
import warnings
from .. import alpaqa as pa
from ..casadi_generator import generate_casadi_control_problem, write_casadi_control_problem_data
from ..cache import get_alpaqa_cache_dir

assert pa.with_casadi_ocp

def _load_casadi_control_problem(sofile, N):
    print("-- Loading:", sofile)
    prob = pa.load_casadi_control_problem(sofile, N=N)
    return prob

def generate_and_compile_casadi_control_problem(
    N: int,
    f: cs.Function,
    l: cs.Function,
    l_N: cs.Function,
    h: cs.Function = None,
    h_N: cs.Function = None,
    c: cs.Function = None,
    c_N: cs.Function = None,
    *,
    U = None,
    D = None,
    D_N = None,
    x_init = None,
    param = None,
    name: str = "alpaqa_control_problem",
    **kwargs,
) -> pa.CasADiControlProblem:
    """Compile the dynamics and cost functions into an alpaqa ControlProblem.

    :param N:    Horizon length.
    :param C:            Bound constraints on u.
    :param D:            Bound constraints on c(x).
    :param D_N:          Bound constraints on c_N(x).
    :param param:        Problem parameter values.
    :param name: Optional string description of the problem (used for filename).
    :param kwargs: Parameters passed to 
                :py:func:`..casadi_generator.generate_casadi_control_problem`.

    :return: Problem specification that can be passed to the solvers.
    """

    cachedir = get_alpaqa_cache_dir()
    cachefile = join(cachedir, 'problems')

    key = base64.b64encode(pickle.dumps(
        (f, l, l_N, h, h_N, c, c_N, name, kwargs))).decode('ascii')

    os.makedirs(cachedir, exist_ok=True)
    with shelve.open(cachefile) as cache:
        if key in cache:
            try:
                uid, soname = cache[key]
                probdir = join(cachedir, str(uid))
                sofile = join(probdir, soname)
                write_casadi_control_problem_data(sofile, U, D, D_N, x_init, param)
                return _load_casadi_control_problem(sofile, N)
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
        cgen = generate_casadi_control_problem(f, l, l_N, h, h_N, c, c_N, name)
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
            warnings.warn(
                f"Multiple compiled CasADi problem files were found for '{name}'"
            )
        sofile = sofile[0]
        soname = os.path.relpath(sofile, probdir)
        cache[key] = uid, soname

        write_casadi_control_problem_data(sofile, U, D, D_N, x_init, param)
        return _load_casadi_control_problem(sofile, N)
