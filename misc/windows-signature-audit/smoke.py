"""Exercise released binaries in fresh processes; retain generated DLLs as evidence."""

import argparse
import json
from pathlib import Path
import subprocess
import sys


CASES = ('symbolic', 'ipopt', 'cvodes', 'jit', 'external')


def run_case(case):
    import casadi as ca

    print(f'CasADi {ca.__version__}: {ca.__file__}', flush=True)
    x = ca.SX.sym('x')
    if case == 'symbolic':
        f = ca.Function('square', [x], [x * x, ca.jacobian(x * x, x)])
        y, dy = f(3)
        assert float(y) == 9 and float(dy) == 6
    elif case == 'ipopt':
        solver = ca.nlpsol('solver', 'ipopt', {'x': x, 'f': (x - 2)**2},
                           {'ipopt.print_level': 0, 'print_time': False})
        result = solver(x0=0)
        assert solver.stats()['success'] and abs(float(result['x']) - 2) < 1e-7
    elif case == 'cvodes':
        import math
        integrator = ca.integrator('integrator', 'cvodes', {'x': x, 'ode': -x}, 0, 1)
        assert abs(float(integrator(x0=1)['xf']) - math.exp(-1)) < 1e-5
    elif case == 'jit':
        f = ca.Function('square', [x], [x * x], {
            'jit': True, 'compiler': 'shell', 'jit_cleanup': False,
            'jit_options': {'cleanup': False, 'verbose': True}})
        assert float(f(3)) == 9
    elif case == 'external':
        f = ca.Function('square', [x], [x * x])
        f.generate('square.c')
        subprocess.run(['cl.exe', '/nologo', '/LD', 'square.c', '/link',
                        '/out:unsigned_probe.dll'], check=True)
        external = ca.external('square', str(Path('unsigned_probe.dll').resolve()))
        assert float(external(3)) == 9
    print(f'PASS: {case}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--case', choices=CASES)
    args = parser.parse_args()
    if args.case:
        run_case(args.case)
        return
    if args.output is None:
        parser.error('--output is required unless --case is specified')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for case in CASES:
        work = output / case
        work.mkdir(exist_ok=True)
        try:
            result = subprocess.run(
                [sys.executable, str(Path(__file__).resolve()), '--case', case],
                cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, errors='replace', timeout=180)
            log, code = result.stdout, result.returncode
        except subprocess.TimeoutExpired as exc:
            log = f'Timed out after 180 seconds\n{exc.stdout!r}\n'
            code = -1
        (work / 'output.txt').write_text(log, encoding='utf-8')
        print(log, flush=True)
        results.append({'case': case, 'exit_code': code})
    (output / 'results.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    sys.exit(any(result['exit_code'] != 0 for result in results))


if __name__ == '__main__':
    main()
