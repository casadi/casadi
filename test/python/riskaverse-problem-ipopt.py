from casadi import SX, vertcat, nlpsol
from casadi.casadi import fmax, fmin
import numpy as np
from numpy import inf

nu = 2
nx = 4
ns = 2
ny = 3
n  = nu + ns + ny
m  = 3

u = SX.sym('u', nu, 1)
s = SX.sym('s', ns, 1)
y = SX.sym('y', ny, 1)
unknowns = vertcat(u, s, y)

Ts = 0.05

A = np.eye(nx)
A[0, 2]           = Ts
A[1, 3]           = Ts
B = np.zeros((nx, nu))
B[2, 0]           = Ts
B[3, 1]           = Ts

x0 = 10 * np.ones((nx, 1))
x1 = A @ x0 + B @ u

Q = 10 * np.eye(nx)
R = 1 * np.eye(nu)

objective = s[0] 
constraints_F1 = [
    y[0] - y[1] - s[0],
    x1.T @ Q @ x1 - s[1],
    x0.T @ Q @ x0 + u.T @ R @ u - (y[0] - y[1] - y[2] - s[1])
]

cs_nlp = {
    'x': unknowns,
    'f': objective,
    'g': vertcat(*constraints_F1),
}
cs_bounds = {
    'lbx': vertcat([-10] * nu + [-inf] * ns + [0] * ny),
    'ubx': vertcat([+10] * nu + [+inf] * ns + [inf] * ny),
    'lbg': vertcat([-inf] * m),
    'ubg': vertcat([0] * m)
}
cs_opts = {
    'verbose': False,
    'ipopt.tol': 1e-5
}

S = nlpsol('S', 'ipopt', cs_nlp, cs_opts)
r = S(x0=np.zeros((n,)), **cs_bounds)
print(f'λ_g = {r["lam_g"]}')
print(f'sol: u = {r["x"][:nu]}')
print(f'     s = {r["x"][nu:nu+ns]}')
print(f'     y = {r["x"][nu+ns:]}')

# Single iteration of PANOC

Σ = np.ones((m,))
y = np.ones((m,))

gxΣ = vertcat(*constraints_F1) + np.divide(y, Σ)
ẑ = fmin(fmax(gxΣ, cs_bounds['lbg']), cs_bounds['ubg'])
proj_sq_sigma_D = (gxΣ - ẑ).T @ np.diag(Σ) @ (gxΣ - ẑ)

ψ = objective + 0.5 * proj_sq_sigma_D

objective = s[0] 
constraints_F1 = [
    y[0] - y[1] - s[0],
    x1.T @ Q @ x1 - s[1],
    x0.T @ Q @ x0 + u.T @ R @ u - (y[0] - y[1] - y[2] - s[1])
]

cs_nlp = {
    'x': unknowns,
    'f': ψ,
    'g': vertcat(),
}
cs_bounds = {
    'lbx': vertcat([-10] * nu + [-inf] * ns + [0] * ny),
    'ubx': vertcat([+10] * nu + [+inf] * ns + [inf] * ny),
    'lbg': vertcat(),
    'ubg': vertcat()
}
cs_opts = {
    'verbose': False,
    'ipopt.tol': 1e-5
}

S = nlpsol('S', 'ipopt', cs_nlp, cs_opts)
r = S(x0=np.zeros((n,)), **cs_bounds)
print(f'λ_g = {r["lam_g"]}')
print(f'sol: u = {r["x"][:nu]}')
print(f'     s = {r["x"][nu:nu+ns]}')
print(f'     y = {r["x"][nu+ns:]}')
