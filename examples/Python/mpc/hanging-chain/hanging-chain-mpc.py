# %% Hanging chain MPC example

import casadi as cs
import numpy as np
from os.path import dirname
import sys

sys.path.append(dirname(__file__))
from hanging_chain_dynamics import HangingChain

# %% Build the model

Ts = 0.05  # Time step [s]
N = 6  # Number of balls
dim = 2  # Dimension (2D or 3D)

model = HangingChain(N, dim, Ts)
y_null, u_null = model.initial_state()  # Initial states and control inputs

model_param = [0.03, 1.6, 0.033 / N]  # Concrete parameters m, D, L

# %% Apply an initial control input to disturb the system

N_dist = 3  # Number of time steps to apply the disturbance for
u_dist = [-0.5, 0.5, 0.5] if dim == 3 else [-0.5, 0.5]  # Disturbance input
y_dist = model.simulate(N_dist, y_null, u_dist, model_param)  # Model states
y_dist = np.hstack((np.array([y_null]).T, y_dist))  # (including initial state)

# %% Simulate the system without a controller

N_sim = 180  # Number of time steps to simulate for
y_sim = model.simulate(N_sim, y_dist[:, -1], u_null, model_param)  # States
y_sim = np.hstack((y_dist, y_sim))  # (including disturbed and initial states)

# %% Define MPC cost and constraints

N_horiz = 12  # MPC horizon length (number of time steps)

y_init = cs.MX.sym("y_init", *y_null.shape)  # Initial state
model_params = cs.MX.sym("params", *model.params.shape)  # Parameters
num_var = dim * N_horiz
U = cs.MX.sym("U", num_var)  # Control signals over horizon
U_mat = model.input_to_matrix(U)  # Input as dim by N_horiz matrix
constr_param = cs.MX.sym("c", 3)  # Coefficients of cubic constraint function
mpc_param = cs.vertcat(y_init, model_params, constr_param)  # All parameters

# Cost

# Stage costs for states and input
stage_y_cost, stage_u_cost = model.generate_cost_funcs()
# Simulate the model with the input over the horizon
mpc_sim = model.simulate(N_horiz, y_init, U_mat, model_params)
# Accumulate the cost of the outputs and inputs
mpc_y_cost = cs.sum2(stage_y_cost.map(N_horiz)(mpc_sim))
mpc_u_cost = cs.sum2(stage_u_cost.map(N_horiz)(U_mat))
mpc_cost = mpc_y_cost + mpc_u_cost

# Constraints

# Cubic constraint function for a single ball in one dimension
g_constr = lambda c, x: c[0] * x**3 + c[1] * x**2 + c[2] * x
# Constraint function for one stage (N balls)
y_c = cs.MX.sym("y_c", y_dist.shape[0])
constr = []
for i in range(N):  # for each ball in the stage except the last,
    yx_n = y_c[dim * i]  # constrain the x, y position of the ball
    yy_n = y_c[dim * i + dim - 1]
    constr += [yy_n - g_constr(constr_param, yx_n)]
constr += [y_c[-1] - g_constr(constr_param, y_c[-dim])]  # Ball N+1
constr_fun = cs.Function("c", [y_c, constr_param], [cs.vertcat(*constr)])
# Constraint function for all stages in the horizon
mpc_constr = cs.vec(constr_fun.map(N_horiz)(mpc_sim, constr_param))
num_constr = (N + 1) * N_horiz
# Fill in the constraint coefficients c(x-a)³ + d(x-a) + b
a, b, c, d = 0.6, -1.4, 5, 2.2
constr_coeff = [c, -3 * a * c, 3 * a * a * c + d]
constr_lb = b - c * a**3 - d * a
# Box constraints on actuator:
C = -1 * np.ones(num_var), +1 * np.ones(num_var)  # lower bound, upper bound
# Constant term of the cubic state constraints as a one-sided box:
D = constr_lb * np.ones(num_constr), +np.inf * np.ones(num_constr)

# Initial parameter value

y_n = np.array(y_dist[:, -1]).ravel()  # Initial state of the chain
n_state = y_n.shape[0]
param_0 = np.concatenate((y_n, model_param, constr_coeff))

# %% NLP formulation

import alpaqa

# Generate C code for the cost and constraint functions, compile them, and load
# them as an alpaqa problem description:
problem = (
    alpaqa.minimize(mpc_cost, U)  # objective and variables         f(x; p)
    .subject_to_box(C)  #           box constraints on variables    x ∊ C
    .subject_to(mpc_constr, D)  #   general constraints             g(x; p) ∊ D
    .with_param(mpc_param)  #       parameter to be changed later   p
    .with_param_value(param_0)  #   initial parameter value
    .with_name(f"hanging_chain_{N_horiz}")  #  name used in generated files
).compile(sym=cs.MX.sym)

# %% NLP solver

from datetime import timedelta

# Configure an alpaqa solver:
solver = alpaqa.ALMSolver(
    alm_params={
        "tolerance": 1e-3,
        "dual_tolerance": 1e-3,
        "initial_penalty": 1e4,
        "max_iter": 100,
        "max_time": timedelta(seconds=0.2),
    },
    inner_solver=alpaqa.PANOCSolver(
        panoc_params={
            "stop_crit": alpaqa.FPRNorm,
            "max_time": timedelta(seconds=0.02),
        },
        lbfgs_params={"memory": N_horiz},
    ),
)

# %% MPC controller


# Wrap the solver in a class that solves the optimal control problem at each
# time step, implementing warm starting:
class MPCController:
    def __init__(self, model: HangingChain, problem: alpaqa.CasADiProblem):
        self.model = model
        self.problem = problem
        self.tot_it = 0
        self.tot_time = timedelta()
        self.max_time = timedelta()
        self.failures = 0
        self.U = np.zeros(problem.num_variables)
        self.λ = np.zeros(problem.num_constraints)

    def __call__(self, y_n):
        y_n = np.array(y_n).ravel()
        # Set the current state as the initial state
        self.problem.param[: y_n.shape[0]] = y_n
        # Shift over the previous solution and Lagrange multipliers
        self.U = np.concatenate((self.U[dim:], self.U[-dim:]))
        self.λ = np.concatenate((self.λ[N + 1 :], self.λ[-N - 1 :]))
        # Solve the optimal control problem
        # (warm start using the shifted previous solution and multipliers)
        self.U, self.λ, stats = solver(self.problem, self.U, self.λ)
        # Print some solver statistics
        print(
            f'{stats["status"]} outer={stats["outer_iterations"]} '
            f'inner={stats["inner"]["iterations"]} time={stats["elapsed_time"]} '
            f'failures={stats["inner_convergence_failures"]}'
        )
        self.tot_it += stats["inner"]["iterations"]
        self.failures += stats["status"] != alpaqa.SolverStatus.Converged
        self.tot_time += stats["elapsed_time"]
        self.max_time = max(self.max_time, stats["elapsed_time"])
        # Print the Lagrange multipliers, shows that constraints are active
        print(np.linalg.norm(self.λ))
        # Return the optimal control signal for the first time step
        return self.model.input_to_matrix(self.U)[:, 0]


# %% Simulate the system using the MPC controller

problem.param = param_0

y_mpc = np.empty((n_state, N_sim))
controller = MPCController(model, problem)
for n in range(N_sim):
    # Solve the optimal control problem:
    u_n = controller(y_n)
    # Apply the first optimal control input to the system and simulate for
    # one time step, then update the state:
    y_n = model.simulate(1, y_n, u_n, model_param).T
    y_mpc[:, n] = y_n
y_mpc = np.hstack((y_dist, y_mpc))

print(f"{controller.tot_it} inner iterations, {controller.failures} failures")
print(
    f"time: {controller.tot_time} (total), {controller.max_time} (max), "
    f"{controller.tot_time / N_sim} (avg)"
)

# %% Visualize the results

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation, patheffects

mpl.rcParams["animation.frame_format"] = "svg"

# Plot the chains
fig, ax = plt.subplots()
x, y, z = model.state_to_pos(y_null)
(line,) = ax.plot(x, y, "-o", label="Without MPC")
(line_ctrl,) = ax.plot(x, y, "-o", label="With MPC")
plt.legend()
plt.ylim([-2.5, 1])
plt.xlim([-0.25, 1.25])

# Plot the state constraints
x = np.linspace(-0.25, 1.25, 256)
y = np.linspace(-2.5, 1, 256)
X, Y = np.meshgrid(x, y)
Z = g_constr(constr_coeff, X) + constr_lb - Y
fx = [patheffects.withTickedStroke(spacing=7, linewidth=0.8)]
cgc = plt.contour(X, Y, Z, [0], colors="tab:green", linewidths=0.8)
cgc.set_path_effects(fx)


class Animation:
    points = []

    def __call__(self, i):
        x, y, z = model.state_to_pos(y_sim[:, i])
        y = z if dim == 3 else y
        for p in self.points:
            p.remove()
        self.points = []
        line.set_xdata(x)
        line.set_ydata(y)
        viol = y - g_constr(constr_coeff, x) + 1e-5 < constr_lb
        if np.sum(viol):
            self.points += ax.plot(x[viol], y[viol], "rx", markersize=12)
        x, y, z = model.state_to_pos(y_mpc[:, i])
        y = z if dim == 3 else y
        line_ctrl.set_xdata(x)
        line_ctrl.set_ydata(y)
        viol = y - g_constr(constr_coeff, x) + 1e-5 < constr_lb
        if np.sum(viol):
            self.points += ax.plot(x[viol], y[viol], "rx", markersize=12)
        return [line, line_ctrl] + self.points


ani = animation.FuncAnimation(
    fig,
    Animation(),
    interval=1000 * Ts,
    blit=True,
    repeat=True,
    frames=1 + N_dist + N_sim,
)

# Show the animation
plt.show()
