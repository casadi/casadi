#%%

import panocpy as pa
import numpy as np
from pprint import pprint
import os
from os.path import join
from tempfile import TemporaryDirectory
import casadi as cs
from pyriskaverse.test.test_optimal_control import TestOptimalControl

rao = TestOptimalControl().constrained_casadi(
    stages=[1, 2], constraint_type="conditional"
)
rao.make_controller()
f = rao.controller.get_function("nlp_f")
g = rao.controller.get_function("nlp_g")
lbx = rao.constraint_bounds["lbx"]
ubx = rao.constraint_bounds["ubx"]
lbg = rao.constraint_bounds["lbg"]
ubg = rao.constraint_bounds["ubg"]

name = "riskmpcproblem"
cgen, n, m, p = pa.generate_casadi_problem(name, f, g)

with TemporaryDirectory(prefix="") as tmpdir:
    cfile = cgen.generate(tmpdir)
    sofile = join(tmpdir, f"{name}.so")
    os.system(f"cc -fPIC -shared -O3 {cfile} -o {sofile}")
    print(sofile)
    prob = pa.load_casadi_problem_with_param(sofile, n, m)

prob.C.lowerbound = lbx
prob.C.upperbound = ubx
prob.D.lowerbound = lbg
prob.D.upperbound = ubg

#%%

panocparams = {
    "max_iter": 2000,
    "print_interval": 200,
    "update_lipschitz_in_linesearch": True,
}

solvers = [
    pa.PANOCSolver(
        pa.PANOCParams(**panocparams),
        pa.LBFGSParams(memory=n),
    ),
    pa.SecondOrderPANOCLBFGSSolver(
        pa.SecondOrderPANOCLBFGSParams(**panocparams),
        pa.LBFGSParams(memory=n),
    ),
]
solverid = 1
almparams = pa.ALMParams(
    max_iter=20,
    print_interval=1,
    preconditioning=False,
    ε_0=1e-2,
    Δ=5,
    Σ_0=1e2,
)
almsolver = pa.ALMSolver(almparams, solvers[solverid])
x0 = np.zeros((n,))
y0 = np.zeros((m,))
p0 = 1e-2 * np.ones((p,))
prob.param = p0
y, x, stats = almsolver(prob, y0, x0)

# print(y)
# print(x)
pprint(stats)

# %%
