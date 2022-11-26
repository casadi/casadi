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

param = [0.03, 1.6, 0.033 / N]  # Concrete parameters m, D, L

# %% Apply an initial control input to disturb the system

N_dist = 3  # Number of time steps to apply the disturbance for
u_dist = [-0.5, 0.5, 0.5] if dim == 3 else [-0.5, 0.5]  # Disturbance input
y_dist = model.simulate(N_dist, y_null, u_dist, param)  # Model states
y_dist = np.hstack((np.array([y_null]).T, y_dist))  # (including initial state)

# %% Simulate the system without a controller

N_sim = 180  # Number of time steps to simulate for
y_sim = model.simulate(N_sim, y_dist[:, -1], u_null, param)  # Model states
y_sim = np.hstack((y_dist, y_sim))  # (including disturbed and initial states)

# %% Define MPC cost and constraints

N_horiz = 12  # MPC horizon length (number of time steps)

y_init = cs.SX.sym("y_init", *y_null.shape)  # Initial state
U = cs.SX.sym("U", dim * N_horiz)  # Control signals over horizon
U_mat = model.input_to_matrix(U)  # Input as dim by N_horiz matrix
constr_param = cs.SX.sym("c", 3)  # Coefficients of cubic constraint function
mpc_param = cs.vertcat(y_init, model.params, constr_param)  # All parameters

# Cost

# Stage costs for states and input
stage_y_cost, stage_u_cost = model.generate_cost_funcs()
# Simulate the model with the input over the horizon
mpc_sim = model.simulate(N_horiz, y_init, U_mat, model.params)
# Accumulate the cost of the outputs and inputs
mpc_y_cost = cs.sum2(stage_y_cost.map(N_horiz)(mpc_sim))
mpc_u_cost = cs.sum2(stage_u_cost.map(N_horiz)(U_mat))
mpc_cost_fun = cs.Function('f_mpc', [U, mpc_param], [mpc_y_cost + mpc_u_cost])

# Constraints

# Cubic constraint function for a single ball in one dimension
g_constr = lambda c, x: c[0] * x**3 + c[1] * x**2 + c[2] * x
# Constraint function for one stage (N balls)
y_c = cs.SX.sym("y_c", y_dist.shape[0])
constr = []
for i in range(N):  # for each ball in the stage except the last,
    yx_n = y_c[dim * i]  # constrain the x, y position of the ball
    yy_n = y_c[dim * i + dim - 1]
    constr += [yy_n - g_constr(constr_param, yx_n)]
constr += [y_c[-1] - g_constr(constr_param, y_c[-dim])]  # Ball N+1
constr_fun = cs.Function("c", [y_c], [cs.vertcat(*constr)])
# Constraint function for all stages in the horizon
mpc_constr = constr_fun.map(N_horiz)(mpc_sim)
mpc_constr_fun = cs.Function("g_mpc", [U, mpc_param], [cs.vec(mpc_constr)])
# Fill in the constraint coefficients c(x-a)³ + d(x-a) + b
a, b, c, d = 0.6, -1.4, 5, 2.2
constr_coeff = [c, -3 * a * c, 3 * a * a * c + d]
constr_lb = b - c * a**3 - d * a

# %% NLP formulation

import alpaqa.casadi_loader as cl

# Generate C code for the cost and constraint function, compile them, and load
# them as an alpaqa problem description:
problem = cl.generate_and_compile_casadi_problem(mpc_cost_fun, mpc_constr_fun)
# Box constraints on actuator:
problem.C.lowerbound = -1 * np.ones((dim * N_horiz, ))
problem.C.upperbound = +1 * np.ones((dim * N_horiz, ))
# Constant term of the cubic state constraints:
problem.D.lowerbound = constr_lb * np.ones((problem.m, ))

# %% NLP solver

import alpaqa as pa
from datetime import timedelta

# Configure an alpaqa solver:
solver = pa.ALMSolver(
    alm_params={
        'ε': 1e-4,
        'δ': 1e-4,
        'Σ_0': 1e4,
        'max_iter': 100,
        'max_time': timedelta(seconds=0.2),
        'max_total_num_retries': 0,
    },
    inner_solver=pa.PANOCSolver(
        panoc_params={
            'stop_crit': pa.ProjGradNorm2,
            'max_time': timedelta(seconds=0.02),
        },
        lbfgs_params={'memory': N_horiz},
    ),
)

# %% MPC controller


# Wrap the solver in a class that solves the optimal control problem at each
# time step, implementing warm starting:
class MPCController:

    def __init__(self, model: HangingChain, problem: pa.CasADiProblem):
        self.model = model
        self.problem = problem
        self.tot_it = 0
        self.tot_time = timedelta()
        self.max_time = timedelta()
        self.failures = 0
        self.U = np.zeros((N_horiz * dim, ))
        self.λ = np.zeros(((N + 1) * N_horiz, ))

    def __call__(self, y_n):
        y_n = np.array(y_n).ravel()
        # Set the current state as the initial state
        self.problem.param[:y_n.shape[0]] = y_n
        # Shift over the previous solution and Lagrange multipliers
        self.U = np.concatenate((self.U[dim:], self.U[-dim:]))
        self.λ = np.concatenate((self.λ[N + 1:], self.λ[-N - 1:]))
        # Solve the optimal control problem
        # (warm start using the shifted previous solution and multipliers)
        self.U, self.λ, stats = solver(self.problem, self.U, self.λ)
        # Print some solver statistics
        print(stats['status'], stats['outer_iterations'],
              stats['inner']['iterations'], stats['elapsed_time'],
              stats['inner_convergence_failures'])
        self.tot_it += stats['inner']['iterations']
        self.failures += stats['status'] != pa.SolverStatus.Converged
        self.tot_time += stats['elapsed_time']
        self.max_time = max(self.max_time, stats['elapsed_time'])
        # Print the Lagrange multipliers, shows that constraints are active
        print(np.linalg.norm(self.λ))
        # Return the optimal control signal for the first time step
        return self.model.input_to_matrix(self.U)[:, 0]


# %% Simulate the system using the MPC controller

y_n = np.array(y_dist[:, -1]).ravel()  # Initial state for controller
n_state = y_n.shape[0]
problem.param = np.concatenate((y_n, param, constr_coeff))

y_mpc = np.empty((n_state, N_sim))
controller = MPCController(model, problem)
for n in range(N_sim):
    # Solve the optimal control problem:
    u_n = controller(y_n)
    # Apply the first optimal control input to the system and simulate for
    # one time step, then update the state:
    y_n = model.simulate(1, y_n, u_n, param).T
    y_mpc[:, n] = y_n
y_mpc = np.hstack((y_dist, y_mpc))

print(f"{controller.tot_it} iterations, {controller.failures} failures")
print(f"time: {controller.tot_time} (total), {controller.max_time} (max), "
      f"{controller.tot_time / N_sim} (avg)")

# %% Visualize the results

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation, patheffects

mpl.rcParams['animation.frame_format'] = 'svg'

# Plot the chains
fig, ax = plt.subplots()
x, y, z = model.state_to_pos(y_null)
line, = ax.plot(x, y, '-o', label='Without MPC')
line_ctrl, = ax.plot(x, y, '-o', label='With MPC')
plt.legend()
plt.ylim([-2.5, 1])
plt.xlim([-0.25, 1.25])

# Plot the state constraints
x = np.linspace(-0.25, 1.25, 256)
y = np.linspace(-2.5, 1, 256)
X, Y = np.meshgrid(x, y)
Z = g_constr(constr_coeff, X) + constr_lb - Y
fx = [patheffects.withTickedStroke(spacing=7, linewidth=0.8)]
cgc = plt.contour(X, Y, Z, [0], colors='tab:green', linewidths=0.8)
plt.setp(cgc.collections, path_effects=fx)


class Animation:
    points = []

    def __call__(self, i):
        x, y, z = model.state_to_pos(y_sim[:, i])
        for p in self.points:
            p.remove()
        self.points = []
        line.set_xdata(x)
        line.set_ydata(y)
        viol = y - g_constr(constr_coeff, x) + 1e-5 < constr_lb
        if np.sum(viol):
            self.points += ax.plot(x[viol], y[viol], 'rx', markersize=12)
        x, y, z = model.state_to_pos(y_mpc[:, i])
        line_ctrl.set_xdata(x)
        line_ctrl.set_ydata(y)
        viol = y - g_constr(constr_coeff, x) + 1e-5 < constr_lb
        if np.sum(viol):
            self.points += ax.plot(x[viol], y[viol], 'rx', markersize=12)
        return [line, line_ctrl] + self.points


ani = animation.FuncAnimation(fig,
                              Animation(),
                              interval=1000 * Ts,
                              blit=True,
                              repeat=True,
                              frames=1 + N_dist + N_sim)

# Show the animation
plt.show()
