import casadi as cs
import numpy as np


def _generate_problem_ss(
    Ts,
    N_hor,
    R_obstacle,
    n_states,
    n_inputs,
    initial_state,
    desired_state,
    f,
    Q,
    R,
    states_lb,
    states_ub,
    inputs_lb,
    inputs_ub,
):
    objective = 0
    state_constr = []
    collision_constr = []

    U = cs.SX.sym("U", n_inputs, N_hor)
    xk = initial_state
    # Forward Euler
    for k in range(N_hor):
        uk = U[:, k]
        err = xk - desired_state
        objective += err.T @ Q @ err + uk.T @ R @ uk
        next_euler = xk + Ts * f(xk, uk)
        if k > 0:
            state_constr = cs.vertcat(state_constr, xk)
            collision_constr = cs.vertcat(
                collision_constr, xk[0:2].T @ xk[0:2]
            )  # ball around (0,0)

        xk = next_euler

    err = xk - desired_state
    objective += 10 * err.T @ Q @ err  # Terminal cost
    state_constr = cs.vertcat(state_constr, xk)
    collision_constr = cs.vertcat(collision_constr, xk[0:2].T @ xk[0:2])

    states_constr_lb = np.tile(states_lb, N_hor)
    states_constr_ub = np.tile(states_ub, N_hor)

    collision_constr_lb = R_obstacle ** 2 * np.ones((N_hor,))
    collision_constr_ub = +np.inf * np.ones((N_hor,))

    nlp = {
        "f": objective,
        "g": cs.vertcat(state_constr, collision_constr),
        "x": U.reshape((N_hor * n_inputs, 1)),
        "p": cs.vertcat(initial_state, desired_state),
    }
    bounds = {
        "lbx": np.tile(inputs_lb, N_hor),
        "ubx": np.tile(inputs_ub, N_hor),
        "lbg": np.concatenate((states_constr_lb, collision_constr_lb)),
        "ubg": np.concatenate((states_constr_ub, collision_constr_ub)),
    }
    first_input_idx = 0

    return f, nlp, bounds, n_states, n_inputs, first_input_idx


def _generate_problem_ms(
    Ts,
    N_hor,
    R_obstacle,
    n_states,
    n_inputs,
    initial_state,
    desired_state,
    f,
    Q,
    R,
    states_lb,
    states_ub,
    inputs_lb,
    inputs_ub,
):
    X, U = cs.SX.sym("X", n_states, N_hor), cs.SX.sym("U", n_inputs, N_hor)
    X0X = cs.horzcat(initial_state, X)
    objective = 0
    dynamics_constr = []
    collision_constr = []

    # Forward Euler
    for k in range(N_hor):
        xk, uk = X0X[:, k], U[:, k]
        xkp1 = X0X[:, k + 1]
        err = xk - desired_state
        objective += err.T @ Q @ err + uk.T @ R @ uk
        next_euler = xk + Ts * f(xk, uk)
        dynamics_constr = cs.vertcat(dynamics_constr, xkp1 - next_euler)
        if k > 0:
            collision_constr = cs.vertcat(
                collision_constr, xk[0:2].T @ xk[0:2]
            )  # ball around (0,0)

    xN = X0X[:, N_hor]
    err = xN - desired_state
    objective += 10 * err.T @ Q @ err  # Terminal cost
    collision_constr = cs.vertcat(collision_constr, xN[0:2].T @ xN[0:2])

    dynamics_constr_lb = np.zeros((N_hor * n_states,))
    dynamics_constr_ub = np.zeros((N_hor * n_states,))

    collision_constr_lb = R_obstacle ** 2 * np.ones((N_hor,))
    collision_constr_ub = +np.inf * np.ones((N_hor,))

    assert dynamics_constr_lb.size == dynamics_constr.size1()
    assert dynamics_constr_ub.size == dynamics_constr.size1()
    assert collision_constr_lb.size == collision_constr.size1()
    assert collision_constr_ub.size == collision_constr.size1()

    nlp = {
        "f": objective,
        "g": cs.vertcat(dynamics_constr, collision_constr),
        "x": cs.vertcat(
            X.reshape((N_hor * n_states, 1)), U.reshape((N_hor * n_inputs, 1))
        ),
        "p": cs.vertcat(initial_state, desired_state),
    }
    bounds = {
        "lbx": np.concatenate((np.tile(states_lb, N_hor), np.tile(inputs_lb, N_hor))),
        "ubx": np.concatenate((np.tile(states_ub, N_hor), np.tile(inputs_ub, N_hor))),
        "lbg": np.concatenate((dynamics_constr_lb, collision_constr_lb)),
        "ubg": np.concatenate((dynamics_constr_ub, collision_constr_ub)),
    }
    first_input_idx = N_hor * n_states

    return f, nlp, bounds, n_states, n_inputs, first_input_idx


# States:
# [0] p_x   longitudinal position
# [1] p_y   lateral position
# [2] ψ     heading angle
# [3] v     longitudinal velocity
#
# Inputs:
# [0] a     longitudinal acceleration
# [1] δ_f   steering angle
def generate_problem(Ts, N_hor, R_obstacle, multipleshooting):

    n_states = 4
    n_inputs = 2

    initial_state = cs.SX.sym("x0", n_states)
    desired_state = cs.SX.sym("xd", n_states)

    # state weights matrix
    Q = cs.diagcat(100, 100, 10, 10)
    # controls weights matrix
    R = cs.diagcat(1e-1, 1e-5)

    # physical constants
    l = 0.5
    l_r = l / 2
    l_f = l / 2

    # state and input constraints
    states_lb = np.array([-10, -10, -np.pi / 3, -3])
    states_ub = np.array([+10, +10, +np.pi / 3, +3])

    inputs_lb = np.array([-5, -np.pi / 4])
    inputs_ub = np.array([+5, +np.pi / 4])

    x, u = cs.SX.sym("x", n_states), cs.SX.sym("u", n_inputs)
    p_x, p_y, ψ, v = cs.vertsplit(x)
    a, δ_f = cs.vertsplit(u)
    β = cs.atan(l_r / (l_f + l_r) * cs.tan(δ_f))
    continuous_dynamics = cs.vertcat(
        v * cs.cos(ψ + β),
        v * cs.sin(ψ + β),
        v / l_r * cs.sin(β),
        a,
    )
    f = cs.Function("f", [x, u], [continuous_dynamics])

    gen = _generate_problem_ms if multipleshooting else _generate_problem_ss
    return gen(
        Ts,
        N_hor,
        R_obstacle,
        n_states,
        n_inputs,
        initial_state,
        desired_state,
        f,
        Q,
        R,
        states_lb,
        states_ub,
        inputs_lb,
        inputs_ub,
    )