## @example mpc/python/bicycle/bicycle-obstacle-avoidance-mpc.py
# This example shows how to call the alpaqa solver from Python code, applied
# to a simple model predictive control problem.

#%% Import necessary libraries and generate the MPC problem

import casadi as cs
import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import timedelta
import os
import sys

sys.path.append(os.path.dirname(__file__))

from bicycle_model import generate_problem

Ts = 0.05  # Sampling time
N_hor = 18  # Horizon length
N_sim = 80 # Length of the simulation

multipleshooting = False
R_obstacle = 2

f, nlp, bounds, n_states, n_inputs, first_input_idx = generate_problem(
    Ts, N_hor, R_obstacle, multipleshooting
)

# %% Build the problem for PANOC+ALM

import panocpy as pa
from tempfile import TemporaryDirectory

name = "mpcproblem"
f_prob = cs.Function("f", [nlp["x"], nlp["p"]], [nlp["f"]])
g_prob = cs.Function("g", [nlp["x"], nlp["p"]], [nlp["g"]])
prob = pa.generate_and_compile_casadi_problem(f_prob, g_prob, name=name)

prob.C.lowerbound = bounds["lbx"]
prob.C.upperbound = bounds["ubx"]
prob.D.lowerbound = bounds["lbg"]
prob.D.upperbound = bounds["ubg"]

#%% PANOC params

lbfgsmem = N_hor
tol = 1e-5
verbose = False

panocparams = {
    "max_iter": 1000,
    "max_time": timedelta(seconds=0.5),
    "print_interval": 10 if verbose else 0,
    "stop_crit": pa.PANOCStopCrit.ProjGradUnitNorm,
    "update_lipschitz_in_linesearch": True,
}

innersolver = pa.PANOCSolver(
    pa.PANOCParams(**panocparams),
    pa.LBFGSParams(memory=lbfgsmem),
)

almparams = pa.ALMParams(
    max_iter=20,
    max_time=timedelta(seconds=1),
    print_interval=1 if verbose else 0,
    preconditioning=False,
    ε=tol,
    δ=tol,
    Δ=5,
    Σ_0=4e5,
    Σ_max=1e12,
)

#%% Simulate MPC

solver = pa.ALMSolver(almparams, innersolver)

state = np.array([-5, 0, 0, 0])
dest = np.array([5, 0.1, 0, 0])

if multipleshooting:
    x_sol = np.concatenate((np.tile(state, N_hor), np.zeros((n_inputs * N_hor,))))
else:
    x_sol = np.zeros((prob.n,))
assert x_sol.size == prob.n
y_sol = np.zeros((prob.m,))


def solve_ocp(state, y_sol, x_sol):
    state = np.reshape(state, (n_states,))
    prob.param = np.concatenate((state, dest))
    t0 = time.perf_counter()
    x_sol, y_sol, stats = solver(problem=prob, x=x_sol, y=y_sol)
    t1 = time.perf_counter()
    return t1 - t0, stats, state, y_sol, x_sol


xs = np.zeros((N_sim, n_states))
times = np.zeros((N_sim,))
for k in range(N_sim):
    t, stats, state, y_sol, x_sol = solve_ocp(state, y_sol, x_sol)
    times[k] = t
    xs[k] = state
    print(
        stats["status"],
        stats["elapsed_time"],
        stats["outer_iterations"],
        stats["inner"]["iterations"],
    )
    input = x_sol[first_input_idx : first_input_idx + n_inputs]
    state += Ts * f(state, input)

#%% Plot

fig_trajectory, ax = plt.subplots(1, 1)
c = plt.Circle((0, 0), R_obstacle)
ax.set_aspect(1)
ax.add_artist(c)
ax.plot(xs[:, 0], xs[:, 1], "r.-")
ax.set_title('Trajectory')

fig_time, ax = plt.subplots(1, 1)
ax.plot(times)
ax.set_title("Run time")
ax.set_xlabel("MPC time step")
ax.set_ylabel("Run time [s]")

plt.show()

# %%
