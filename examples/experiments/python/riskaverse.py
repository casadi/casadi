#%%

import panocpy as pa
import numpy as np
from pprint import pprint
import os
from os.path import join
from tempfile import TemporaryDirectory
import casadi as cs
from pyriskaverse.test.test_optimal_control import TestOptimalControl
import opengen as og

np.random.seed(0)

rao = TestOptimalControl().constrained_casadi(
    stages=[1, 2], constraint_type="conditional"
)
# rao.aux_vars_formulation_cost = rao.ONLY_S
rao.make_controller()
f = rao.controller.get_function("nlp_f")
g = rao.controller.get_function("nlp_g")
lbx = rao.constraint_bounds["lbx"]
ubx = rao.constraint_bounds["ubx"]
lbg = rao.constraint_bounds["lbg"]
ubg = rao.constraint_bounds["ubg"]
x_sym = f.sx_in(0)
p_sym = f.sx_in(1)

name = "riskmpcproblem"
cgen, n, m, num_p = pa.generate_casadi_problem(name, f, g)

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
    "print_interval": 10,
    "update_lipschitz_in_linesearch": False,
    "stop_crit": pa.PANOCStopCrit.ProjGradNorm,
}
lbfgsmem = 10

solvers = [
    pa.PANOCSolver(
        pa.PANOCParams(**panocparams),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
    pa.StructuredPANOCLBFGSSolver(
        pa.StructuredPANOCLBFGSParams(**panocparams, nonmonotone_linesearch=0),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
]
solverid = 1
almparams = pa.ALMParams(
    max_iter=20,
    print_interval=1,
    preconditioning=False,
    ε=1e-7,
    δ=1e-4,
    ε_0=1e-2,
    Δ=5,
    Σ_0=0,
)
almsolver = pa.ALMSolver(almparams, solvers[solverid])
x0 = np.zeros((n,))
y0 = np.zeros((m,))
p0 = 1e-3 * np.ones((num_p,))
prob.param = p0
y, x, stats = almsolver(prob, y0, x0)

# print(y)
# print(x)
print(f"inner iterations: {stats['inner']['iterations']}")
print(
    f"time: {stats['elapsed_time'].seconds * 1e3 + stats['elapsed_time'].microseconds / 1000} ms"
)
print(f(x, p0))

# %%

sol = rao.control(p0)
print(sol.status_code)
print(sol.states)
print(sol.control_actions)
print(sol.optimal_value)
print(sol.status_msg)
x_ipopt = np.array(sol.raw_output["x"])

#%% OpEn

C = og.constraints.Rectangle(xmin=[*np.array(lbx)], xmax=[*np.array(ubx)])
D = og.constraints.Rectangle(xmin=[*np.array(lbg)], xmax=[*np.array(ubg)])

problem = (
    og.builder.Problem(x_sym, p_sym, f(x_sym, p_sym))
    .with_constraints(C)
    .with_aug_lagrangian_constraints(g(x_sym, p_sym), D)
)

meta = (
    og.config.OptimizerMeta()
    .with_version("0.0.0")
    .with_authors(["Pieter Pas"])
    .with_optimizer_name(name)
)
build_config = (
    og.config.BuildConfiguration()
    .with_build_directory("/tmp/open")
    .with_build_mode("release")
    .with_tcp_interface_config()
)
solver_config = (
    og.config.SolverConfiguration()
    .with_lbfgs_memory(lbfgsmem)
    .with_tolerance(almparams.ε)
    .with_delta_tolerance(almparams.δ)
    .with_max_inner_iterations(10000)
    .with_penalty_weight_update_factor(almparams.Δ)
    # .with_initial_penalty(almparams.Σ_0)
    .with_initial_tolerance(almparams.ε_0)
    .with_sufficient_decrease_coefficient(almparams.θ)
    .with_inner_tolerance_update_factor(almparams.ρ)
    .with_max_outer_iterations(almparams.max_iter)
    .with_max_duration_micros(
        almparams.max_time.microseconds + 1e6 * almparams.max_time.seconds
    )
)
builder = og.builder.OpEnOptimizerBuilder(
    problem,
    metadata=meta,
    build_configuration=build_config,
    solver_configuration=solver_config,
)
builder.build()

#%%

mng = og.tcp.OptimizerTcpManager(join("/tmp/open", name))
mng.start()

pong = mng.ping()  # check if the server is alive
print(pong)
response = mng.call(p0)  # call the solver over TCP

if response.is_ok():
    # Solver returned a solution
    solution_data = response.get()
    u_star = solution_data.solution
    exit_status = solution_data.exit_status
    solver_time = solution_data.solve_time_ms
    pprint(solution_data.__dict__)
else:
    # Invocation failed - an error report is returned
    solver_error = response.get()
    print(f"{solver_error.code}: {solver_error.message}")

mng.kill()

# %%

import matplotlib.pyplot as plt

plt.figure()
plt.semilogy(np.abs(x_ipopt))
plt.semilogy(np.abs(x))

plt.figure()
plt.plot(x_ipopt)
plt.plot(x)
plt.show()
# %%
