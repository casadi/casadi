from typing import Union
import casadi as cs
import numpy as np
from casadi import vertcat as vc


class HangingChain:

    def __init__(self, N: int, dim: int, Ts: float = 0.05):
        self.N = N
        self.dim = dim

        self.y1 = cs.SX.sym("y1", dim * N, 1)  # state: balls 1→N positions
        self.y2 = cs.SX.sym("y2", dim * N, 1)  # state: balls 1→N velocities
        self.y3 = cs.SX.sym("y3", dim, 1)  # state: ball  1+N position
        self.u = cs.SX.sym("u", dim, 1)  # input: ball  1+N velocity
        self.y = vc(self.y1, self.y2, self.y3)  # full state vector

        self.m = cs.SX.sym("m")  # mass
        self.D = cs.SX.sym("D")  # spring constant
        self.L = cs.SX.sym("L")  # spring length
        self.params = vc(self.m, self.D, self.L)

        self.g = np.array([0, 0, -9.81] if dim == 3 else [0, -9.81])  # gravity
        self.x0 = np.zeros((dim, ))  # ball 0 position
        self.x_end = np.eye(1, dim, 0).ravel()  # ball N+1 reference position

        self._build_dynamics(Ts)

    def _build_dynamics(self, Ts):
        y, y1, y2, y3, u = self.y, self.y1, self.y2, self.y3, self.u
        dist = lambda xa, xb: cs.norm_2(xa - xb)
        N, d = self.N, self.dim
        p = self.params

        # Continuous-time dynamics y' = f(y, u; p)

        f1 = [y2]
        f2 = []
        for i in range(N):
            xi = y1[d * i:d * i + d]
            xip1 = y1[d * i + d:d * i + d * 2] if i < N - 1 else y3
            Fiip1 = self.D * (1 - self.L / dist(xip1, xi)) * (xip1 - xi)
            xim1 = y1[d * i - d:d * i] if i > 0 else self.x0
            Fim1i = self.D * (1 - self.L / dist(xi, xim1)) * (xi - xim1)
            fi = (Fiip1 - Fim1i) / self.m + self.g
            f2 += [fi]
        f3 = [u]

        f_expr = vc(*f1, *f2, *f3)
        self.f = cs.Function("f", [y, u, p], [f_expr], ["y", "u", "p"], ["y'"])

        # Discretize dynamics y[k+1] = f_d(y[k], u[k]; p)

        # 4th order Runge-Kutta integrator
        k1 = self.f(y, u, p)
        k2 = self.f(y + Ts * k1 / 2, u, p)
        k3 = self.f(y + Ts * k2 / 2, u, p)
        k4 = self.f(y + Ts * k3, u, p)

        # Discrete-time dynamics
        f_d_expr = y + (Ts / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        self.f_d = cs.Function("f", [y, u, p], [f_d_expr])

        return self.f_d

    def state_to_pos(self, y):
        N, d = self.N, self.dim
        rav = lambda x: np.array(x).ravel()
        xdim = lambda y, i: np.concatenate(
            ([0], rav(y[i:d * N:d]), rav(y[-d + i])))
        if d == 3:
            return (xdim(y, 0), xdim(y, 1), xdim(y, 2))
        else:
            return (xdim(y, 0), xdim(y, 1), np.zeros((N + 1, )))

    def input_to_matrix(self, u):
        """
        Reshape the input signal from a vector into a dim × N_horiz matrix (note
        that CasADi matrices are stored column-wise and NumPy arrays row-wise)
        """
        if isinstance(u, np.ndarray):
            return u.reshape((self.dim, u.shape[0] // self.dim), order='F')
        else:
            return u.reshape((self.dim, u.shape[0] // self.dim))

    def simulate(self, N_sim: int, y_0: np.ndarray,
                 u: Union[np.ndarray, list, cs.SX.sym, cs.MX.sym],
                 p: Union[np.ndarray, list, cs.SX.sym, cs.MX.sym]):
        if isinstance(u, list):
            u = np.array(u)
        if isinstance(u, np.ndarray):
            if u.ndim == 1 or (u.ndim == 2 and u.shape[1] == 1):
                if u.shape[0] == self.dim:
                    u = np.tile(u, (N_sim, 1)).T
        return self.f_d.mapaccum(N_sim)(y_0, u, p)

    def initial_state(self):
        N, d = self.N, self.dim
        y1_0 = np.zeros((d * N))
        y1_0[0::d] = np.arange(1, N + 1) / (N + 1)
        y2_0 = np.zeros((d * N))
        y3_0 = np.zeros((d, ))
        y3_0[0] = 1

        y_null = np.concatenate((y1_0, y2_0, y3_0))
        u_null = np.zeros((d, ))

        return y_null, u_null

    def generate_cost_funcs(self, α=25, β=1, γ=0.01):
        N, d = self.N, self.dim
        y1t = cs.SX.sym("y1t", d * N, 1)
        y2t = cs.SX.sym("y2t", d * N, 1)
        y3t = cs.SX.sym("y3t", d, 1)
        ut = cs.SX.sym("ut", d, 1)
        yt = cs.vertcat(y1t, y2t, y3t)

        L_cost_x = α * cs.sumsqr(y3t - self.x_end)
        for i in range(N):
            xdi = y2t[d * i:d * i + d]
            L_cost_x += β * cs.sumsqr(xdi)
        L_cost_u = γ * cs.sumsqr(ut)
        return cs.Function("L_cost_x", [yt], [L_cost_x]), \
            cs.Function("L_cost_u", [ut], [L_cost_u])
