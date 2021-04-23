# Based on https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi/blob/master/workshop_github/Python_Implementation/mpc_code.py

import casadi as cs
import numpy as np
from casadi import sin, cos, pi, Function, CodeGenerator, gradient, jtimes, hessian
from sys import argv

if len(argv) < 2:
    print(f"Usage:    {argv[0]} <name>")
    exit(1)

# setting matrix_weights' variables
Q_x = 100
Q_y = 100
Q_theta = 2000
R1 = 1
R2 = 1
R3 = 1
R4 = 1

step_horizon = 0.1  # time between steps in seconds
N = 10              # number of look ahead steps
rob_diam = 0.3      # diameter of the robot
wheel_radius = 1    # wheel radius
Lx = 0.3            # L in J Matrix (half robot x-axis length)
Ly = 0.3            # l in J Matrix (half robot y-axis length)
sim_time = 200      # simulation time

# specs
x_init = 0
y_init = 0
theta_init = 0
x_target = 15
y_target = 10
theta_target = pi/4

v_max = 1
v_min = -1


def shift_timestep(step_horizon, t0, state_init, u, f):
    f_value = f(state_init, u[:, 0])
    next_state = cs.DM.full(state_init + (step_horizon * f_value))

    t0 = t0 + step_horizon
    u0 = cs.horzcat(
        u[:, 1:],
        cs.reshape(u[:, -1], -1, 1)
    )

    return t0, next_state, u0


def DM2Arr(dm):
    return np.array(dm.full())


# state symbolic variables
x = cs.SX.sym('x')
y = cs.SX.sym('y')
theta = cs.SX.sym('theta')
states = cs.vertcat(
    x,
    y,
    theta
)
n_states = states.numel()

# control symbolic variables
V_a = cs.SX.sym('V_a')
V_b = cs.SX.sym('V_b')
V_c = cs.SX.sym('V_c')
V_d = cs.SX.sym('V_d')
controls = cs.vertcat(
    V_a,
    V_b,
    V_c,
    V_d
)
n_controls = controls.numel()

# matrix containing all states over all time steps +1 (each column is a state vector)
X = cs.SX.sym('X', n_states, N + 1)

# matrix containing all control actions over all time steps (each column is an action vector)
U = cs.SX.sym('U', n_controls, N)

# column vector for storing initial state and target state
# P = cs.SX.sym('P', n_states + n_states)
x_init = 0
y_init = 0
theta_init = 0

x_target = 15
y_target = 10
theta_target = pi/4

P =  cs.vertcat(x_init, 
                y_init,
                theta_init,
                x_target,
                y_target,
                theta_target)

# state weights matrix (Q_X, Q_Y, Q_THETA)
Q = cs.diagcat(Q_x, Q_y, Q_theta)

# controls weights matrix
R = cs.diagcat(R1, R2, R3, R4)

# discretization model (e.g. x2 = f(x1, v, t) = x1 + v * dt)
rot_3d_z = cs.vertcat(
    cs.horzcat(cos(theta), -sin(theta), 0),
    cs.horzcat(sin(theta),  cos(theta), 0),
    cs.horzcat(         0,           0, 1)
)
# Mecanum wheel transfer function which can be found here: 
# https://www.researchgate.net/publication/334319114_Model_Predictive_Control_for_a_Mecanum-wheeled_robot_in_Dynamical_Environments
J = (wheel_radius/4) * cs.DM([
    [         1,         1,          1,         1],
    [        -1,         1,          1,        -1],
    [-1/(Lx+Ly), 1/(Lx+Ly), -1/(Lx+Ly), 1/(Lx+Ly)]
])
# RHS = states + J @ controls * step_horizon  # Euler discretization
RHS = rot_3d_z @ J @ controls
# maps controls from [va, vb, vc, vd].T to [vx, vy, omega].T
f = cs.Function('f', [states, controls], [RHS])


cost_fn = 0  # cost function
g = X[:, 0] - P[:n_states]  # constraints in the equation


# runge kutta
for k in range(N):
    st = X[:, k]
    con = U[:, k]
    cost_fn = cost_fn \
        + (st - P[n_states:]).T @ Q @ (st - P[n_states:]) \
        + con.T @ R @ con
    st_next = X[:, k+1]
    k1 = f(st, con)
    k2 = f(st + step_horizon/2*k1, con)
    k3 = f(st + step_horizon/2*k2, con)
    k4 = f(st + step_horizon * k3, con)
    st_next_RK4 = st + (step_horizon / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    g = cs.vertcat(g, st_next - st_next_RK4)


OPT_variables = cs.vertcat(
    X.reshape((-1, 1)),   # Example: 3x11 ---> 33x1 where 3=states, 11=N+1
    U.reshape((-1, 1))
)

lbx = cs.DM.zeros((n_states*(N+1) + n_controls*N, 1))
ubx = cs.DM.zeros((n_states*(N+1) + n_controls*N, 1))

lbx[0: n_states*(N+1): n_states] = -cs.inf     # X lower bound
lbx[1: n_states*(N+1): n_states] = -cs.inf     # Y lower bound
lbx[2: n_states*(N+1): n_states] = -cs.inf     # theta lower bound

ubx[0: n_states*(N+1): n_states] = cs.inf      # X upper bound
ubx[1: n_states*(N+1): n_states] = cs.inf      # Y upper bound
ubx[2: n_states*(N+1): n_states] = cs.inf      # theta upper bound

lbx[n_states*(N+1):] = v_min                  # v lower bound for all V
ubx[n_states*(N+1):] = v_max                  # v upper bound for all V

lbg = cs.DM.zeros((n_states*(N+1), 1))  # constraints lower bound
ubg = cs.DM.zeros((n_states*(N+1), 1))  # constraints upper bound

args = {
    'lbg': lbg,
    'ubg': ubg,
    'lbx': lbx,
    'ubx': ubx
}

t0 = 0
state_init = cs.DM([x_init, y_init, theta_init])        # initial state
state_target = cs.DM([x_target, y_target, theta_target])  # target state

# xx = DM(state_init)
t = cs.DM(t0)

u0 = cs.DM.zeros((n_controls, N))  # initial control
X0 = cs.repmat(state_init, 1, N+1)         # initial state full


mpc_iter = 0
cat_states = DM2Arr(X0)
cat_controls = DM2Arr(u0[:, 0])
times = np.array([[0]])


###############################################################################

unknwns = OPT_variables
n_unknwns = unknwns.shape[0]
m_constr = g.shape[0]

w = cs.SX.sym("w", m_constr)
λ = cs.SX.sym("λ", m_constr)
v = cs.SX.sym("v", n_unknwns)
f = cost_fn

L = f + λ.T @ g

cg = CodeGenerator(f"{argv[1]}.c")
cg.add(Function("f", [unknwns],
                [f],
                ["x"], ["f"]))
cg.add(Function("grad_f", [unknwns],
                [gradient(f, unknwns)],
                ["x"], ["grad_f"]))
cg.add(Function("g", [unknwns],
                [g],
                ["x"], ["g"]))
cg.add(Function("grad_g", [unknwns, w],
                [jtimes(g, unknwns, w, True)],
                ["x", "w"], ["grad_g"]))
cg.add(Function("hess_L", [unknwns, λ],
                [hessian(L, unknwns)[0]],
                ["x", "y"], ["hess_L"]))
cg.add(Function("hess_L_prod", [unknwns, λ, v],
                [gradient(jtimes(L, unknwns, v, False), unknwns)],
                ["x", "y", "v"], ["hess_L_prod"]))
cg.generate()

solver_x0 = cs.vertcat(
    cs.reshape(X0, n_states*(N+1), 1),
    cs.reshape(u0, n_controls*N, 1)
)
np.savetxt(f'{argv[1]}-lbx.csv', np.array(lbx))
np.savetxt(f'{argv[1]}-ubx.csv', np.array(ubx))
np.savetxt(f'{argv[1]}-lbg.csv', np.array(lbg))
np.savetxt(f'{argv[1]}-ubg.csv', np.array(ubg))