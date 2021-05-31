#%%
import casadi as cs
import numpy as np
import time
import matplotlib.pyplot as plt
from bicycle_model import generate_problem

Ts = 0.05  # Sampling time
N_hor = 16  # Horizon length
N_sim = 60

multipleshooting = False
tol = 1e-4
R_obstacle = 2

f, nlp, bounds, n_states, n_inputs, first_input_idx = generate_problem(Ts, N_hor, R_obstacle, multipleshooting)

opts = {
    "verbose_init": False,
    "print_time": False,
    "ipopt.tol": tol,
    "ipopt.print_level": 0,
}

#%% NLPSOL

solver = cs.nlpsol("solver", "ipopt", nlp, opts)
p0 = np.array([*[-5, 0, 0, 0], *[0, 0, 0, 0]])  # initial state, desired state
sol = solver(p=p0, **bounds)

#%% Simulate

xs = np.zeros((N_sim, n_states))
ipopt_times = np.zeros((N_sim,))

state = np.array([-5, 0, 0, 0])
dest = np.array([5, 0.1, 0, 0])
for k in range(N_sim):
    state = np.reshape(state, (n_states,))
    xs[k] = state
    p = np.concatenate((state, dest))
    t0 = time.perf_counter()
    sol = solver(p=p, **bounds)
    t1 = time.perf_counter()
    ipopt_times[k] = t1 - t0
    input = sol["x"][first_input_idx : first_input_idx + 2]
    state += Ts * f(state, input)

#%% Plot

fig, (ax_ipopt, ax_panoc, ax_open) = plt.subplots(3, 1)
c = plt.Circle((0, 0), R_obstacle)
ax_ipopt.set_aspect(1)
ax_ipopt.add_artist(c)
ax_ipopt.plot(xs[:, 0], xs[:, 1], "r.-")

# %%

import panocpy as pa
from tempfile import TemporaryDirectory
import os
from os.path import join

name = "mpcproblem"
f_prob = cs.Function("f", [nlp["x"], nlp["p"]], [nlp["f"]])
g_prob = cs.Function("g", [nlp["x"], nlp["p"]], [nlp["g"]])
cgen, n, m, num_p = pa.generate_casadi_problem(name, f_prob, g_prob)

with TemporaryDirectory(prefix="") as tmpdir:
    cfile = cgen.generate(tmpdir)
    sofile = join(tmpdir, f"{name}.so")
    os.system(f"cc -fPIC -shared -O3 -march=native {cfile} -o {sofile}")
    print(sofile)
    prob = pa.load_casadi_problem_with_param(sofile, n, m)

prob.C.lowerbound = bounds["lbx"]
prob.C.upperbound = bounds["ubx"]
prob.D.lowerbound = bounds["lbg"]
prob.D.upperbound = bounds["ubg"]

#%% PANOC params
verbose = False
from datetime import timedelta
panocparams = {
    "max_iter": 1000,
    "print_interval": 1000 if verbose else 0,
    "stop_crit": pa.PANOCStopCrit.ApproxKKT,
    "update_lipschitz_in_linesearch": False,
}
lbfgsmem = N_hor * 4

solvers = [
    pa.PANOCSolver(
        pa.PANOCParams(**panocparams),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
    pa.SecondOrderPANOCLBFGSSolver(
        pa.SecondOrderPANOCLBFGSParams(**panocparams),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
]
solverid = 1
almparams = pa.ALMParams(
    max_iter=100,
    max_time=timedelta(seconds=10),
    print_interval=1 if verbose else 0,
    preconditioning=False,
    ε=tol,
    δ=tol,
    ε_0=1e0,
    Δ=10,
    Σ_0=(4e5 if multipleshooting else 4e5),
    max_total_num_retries=0,
    σ_0=1e-1,
    Σ_max=1e12,
)
almsolver = pa.ALMSolver(almparams, solvers[solverid])
y0 = np.zeros((m,))

#%% Print fpr

# N_sim = 40
xs = np.zeros((N_sim, n_states))
panoc_times = np.zeros((N_sim,))
avg_τs = np.zeros((N_sim,))
panoc_inner_iters = np.zeros((N_sim,))

state = np.array([-5, 0, 0, 0])
dest = np.array([5, 0.1, 0, 0])
if multipleshooting:
    x_sol = np.concatenate((np.tile(state, N_hor), np.zeros((n_inputs * N_hor,))))
else:
    x_sol = np.zeros((n,))
assert x_sol.size == n
y_sol = np.zeros((m,))

for k in range(0):
    state = np.reshape(state, (n_states,))
    xs[k] = state
    prob.param = np.concatenate((state, dest))
    t0 = time.perf_counter()
    y_sol, x_sol, stats = almsolver(prob, y_sol, x_sol)
    t1 = time.perf_counter()
    panoc_times[k] = t1 - t0
    # pprint(stats)
    avg_τs[k] = stats['inner']['sum_τ'] / stats['inner']['count_τ'] if stats['inner']['count_τ'] > 0 else 1
    panoc_inner_iters[k] = stats['inner']['iterations']
    # print(stats["status"], stats["elapsed_time"], stats['outer_iterations'], stats['inner']['iterations'], avg_τs[k])
    input = x_sol[first_input_idx : first_input_idx + 2]
    state += Ts * f(state, input)

class fpr_logger:
    def __init__(self):
        self.alm_it = -1
        self.data = []
    
    def update(self, s):
        if s.k == 0:
            self.alm_it += 1
            self.data.append([])
        self.data[self.alm_it].append(s.fpr)

logger = fpr_logger()

state = np.reshape(state, (n_states,))
prob.param = np.concatenate((state, dest))
almsolver.inner_solver().set_progress_callback(logger.update)
y_sol, x_sol, stats = almsolver(prob, y_sol, x_sol)

from pprint import pprint

pprint(logger.data)
plt.figure()
for i, fprs in enumerate(logger.data):
    plt.semilogy(fprs, '.-', label=str(i))
plt.legend()
plt.savefig(f'fpr-upd-ls-cond={panocparams["update_lipschitz_in_linesearch"]}.pdf')

#%% Simulate

from pprint import pprint

# N_sim = 40
xs = np.zeros((N_sim, n_states))
panoc_times = np.zeros((N_sim,))
avg_τs = np.zeros((N_sim,))
panoc_inner_iters = np.zeros((N_sim,))

state = np.array([-5, 0, 0, 0])
dest = np.array([5, 0.1, 0, 0])
if multipleshooting:
    x_sol = np.concatenate((np.tile(state, N_hor), np.zeros((n_inputs * N_hor,))))
else:
    x_sol = np.zeros((n,))
assert x_sol.size == n
y_sol = np.zeros((m,))

for k in range(N_sim):
    state = np.reshape(state, (n_states,))
    xs[k] = state
    prob.param = np.concatenate((state, dest))
    t0 = time.perf_counter()
    y_sol, x_sol, stats = almsolver(prob, y_sol, x_sol)
    t1 = time.perf_counter()
    panoc_times[k] = t1 - t0
    # pprint(stats)
    avg_τs[k] = stats['inner']['sum_τ'] / stats['inner']['count_τ'] if stats['inner']['count_τ'] > 0 else 1
    panoc_inner_iters[k] = stats['inner']['iterations']
    print(stats["status"], stats["elapsed_time"], stats['outer_iterations'], stats['inner']['iterations'], avg_τs[k])
    input = x_sol[first_input_idx : first_input_idx + 2]
    state += Ts * f(state, input)

print(panoc_inner_iters.sum())

#%% Plot

import matplotlib.pyplot as plt

c = plt.Circle((0, 0), R_obstacle)
ax_panoc.set_aspect(1)
ax_panoc.add_artist(c)
ax_panoc.plot(xs[:, 0], xs[:, 1], "r.-")
# %%

_, (ax_1, ax_2) = plt.subplots(2, 1)
ax_1.plot(avg_τs, ".-", label="avg τ")
ax_1.set_ylabel('avg τ')
ax_2.plot(panoc_inner_iters, ".-", label="iter")
ax_2.set_ylabel('Inner iter')
ax_2.set_xlabel('Time step')

# %%
import opengen as og

C = og.constraints.Rectangle(xmin=[*prob.C.lowerbound], xmax=[*prob.C.upperbound])
D = og.constraints.Rectangle(xmin=[*prob.D.lowerbound], xmax=[*prob.D.upperbound])

problem = (
    og.builder.Problem(nlp['x'], nlp['p'], nlp['f'])
    .with_constraints(C)
    .with_aug_lagrangian_constraints(nlp['g'], D)
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
    .with_max_inner_iterations(almparams.max_iter * panocparams['max_iter'])
    .with_penalty_weight_update_factor(almparams.Δ)
    .with_initial_penalty(almparams.Σ_0)
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

xs = np.zeros((N_sim, n_states))
open_times = np.zeros((N_sim,))
open_inner_iters = np.zeros((N_sim,))

state = np.array([-5, 0, 0, 0])
dest = np.array([5, 0.1, 0, 0])
if multipleshooting:
    x_sol = np.concatenate((np.tile(state, N_hor), np.zeros((n_inputs * N_hor,))))
else:
    x_sol = np.zeros((n,))
assert x_sol.size == n
y_sol = np.zeros((m,))

for k in range(N_sim):
    state = np.reshape(state, (n_states,))
    xs[k] = state
    t0 = time.perf_counter()
    stats = mng.call(np.concatenate((state, dest)), initial_guess=x_sol, initial_y=y_sol)
    t1 = time.perf_counter()
    open_times[k] = t1 - t0

    if stats.is_ok():
        # Solver returned a solution
        solution_data = stats.get()
        x_sol = solution_data.solution
        y_sol = solution_data.lagrange_multipliers
        print(solution_data.exit_status, solution_data.solve_time_ms, solution_data.num_outer_iterations, solution_data.num_inner_iterations)
        open_inner_iters[k] = solution_data.num_inner_iterations
    input = x_sol[first_input_idx : first_input_idx + 2]
    state += Ts * f(state, input)

mng.kill()

# %%

c = plt.Circle((0, 0), R_obstacle)
ax_open.set_aspect(1)
ax_open.add_artist(c)
ax_open.plot(xs[:, 0], xs[:, 1], "r.-")

# %%

plt.figure()
plt.semilogy(panoc_inner_iters, ".-", label="PANOC+ALM")
plt.semilogy(open_inner_iters, ".-", label="OpEn")
plt.title('Inner iterations')
plt.xlabel('Time step')
plt.legend()

# %%

plt.figure()
plt.semilogy(ipopt_times, ".-", label="Ipopt")
plt.semilogy(panoc_times, ".-", label="PANOC+ALM")
plt.semilogy(open_times, ".-", label="OpEn")
plt.title('Time')
plt.xlabel('Time step')
plt.legend()

# %%
plt.show()

# %%



# %%
