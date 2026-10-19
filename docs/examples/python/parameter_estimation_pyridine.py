#
#     MIT No Attribution
#
#     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# Rate constant estimation for the denitrogenation of pyridine, after H.G. Bock,
# "Numerical treatment of inverse problems in chemical reaction kinetics" (1981):
# multiple shooting and a generalized Gauss-Newton method, compared to single shooting.
# The 'rosenbrock' integrator returns the exact derivatives of its adaptive discretization,
# so the Gauss-Newton iterations behave as predicted by the theory.
#
import casadi as ca
import numpy as np

# Reaction kinetics: pyridine, piperidine, pentylamine, N-pentylpiperidine, dipentylamine,
# ammonia, pentane
c = ca.MX.sym('c', 7)
p = ca.MX.sym('p', 11)
A, B, C, D, E, F, G = ca.vertsplit(c)
rhs = ca.vertcat(
  -p[0]*A + p[8]*B,
  p[0]*A - p[1]*B - p[2]*C*B + p[6]*D - p[8]*B + p[9]*D*F,
  p[1]*B - p[2]*B*C - 2*p[3]*C*C - p[5]*C + p[7]*E + p[9]*D*F + 2*p[10]*E*F,
  p[2]*B*C - p[4]*D - p[6]*D - p[9]*D*F,
  p[3]*C*C + p[4]*D - p[7]*E - p[10]*E*F,
  p[2]*B*C + p[3]*C*C + p[5]*C - p[9]*D*F - p[10]*E*F,
  p[5]*C + p[6]*D + p[7]*E)
dae = {'x': c, 'p': p, 'ode': rhs}
opts = {'reltol': 1e-8, 'abstol': 1e-10}

# Synthetic measurements of all species at t = 0, 0.5, ..., 5.5, using the rate constants
# estimated in the paper
p_true = [1.81, 0.894, 29.4, 9.21, 0.0580, 2.43, 0.0644, 5.55, 0.0201, 0.577, 2.15]
c0_true = [1, 0, 0, 0, 0, 0, 0]
t_meas = np.arange(12) * 0.5
sim = ca.integrator('sim', 'rosenbrock', dae, 0, t_meas[1:], opts)
c_true = np.hstack([np.array(c0_true)[:, None], np.array(sim(x0=c0_true, p=p_true)['xf'])])
sigma = 0.001
data = c_true + sigma * np.random.default_rng(1).standard_normal(c_true.shape)


class Linearization:
  """Generalized inverse J^+ of the linearized problem min |F1 + J1 dx|^2, F2 + J2 dx = 0,
  by the null-space method: J2^T = [Q1 Q2] [R; 0], dx = Q1 a + Q2 b"""
  def __init__(self, J1, J2):
    m = J2.shape[0]
    Q, R = np.linalg.qr(J2.T, mode='complete')
    self.Q1, self.Q2, self.R = Q[:, :m], Q[:, m:], R[:m, :]
    self.J1, self.A = J1, J1 @ self.Q2

  def step(self, f1, f2):
    a = -np.linalg.solve(self.R.T, f2) if f2.size else np.zeros(0)
    b = np.linalg.lstsq(self.A, -(f1 + self.J1 @ self.Q1 @ a), rcond=None)[0]
    return self.Q1 @ a + self.Q2 @ b

  def covariance(self):
    # Covariance of the estimate for unit-variance residuals, cf. (8.41)
    Ap = np.linalg.pinv(self.A)
    return self.Q2 @ Ap @ Ap.T @ self.Q2.T


def gauss_newton(F1, F2, x0, label, max_iter=50, tol=1e-10):
  """Generalized Gauss-Newton for min |F1(x)|^2 s.t. F2(x) = 0, damped with the natural
  level function |J_k^+ F(x)| (restricted monotonicity test)"""
  x = ca.MX.sym('x', x0.size)
  r1, r2 = F1(x), F2(x)
  fun = ca.Function('fun', [x], [r1, ca.jacobian(r1, x), r2, ca.jacobian(r2, x)])
  res = ca.Function('res', [x], [r1, r2])
  def residuals(xk):
    try:
      return [np.array(e).ravel() for e in res(xk)]
    except RuntimeError:
      return None  # integration failed, e.g. for wild rate constants
  xk = x0.copy()
  lam = 1.
  for it in range(max_iter):
    f1, J1, f2, J2 = [np.array(e) for e in fun(xk)]
    f1, f2 = f1.ravel(), f2.ravel()
    lin = Linearization(J1, J2)
    dx = lin.step(f1, f2)
    # Damping: accept lam if the simplified increment at the trial point is smaller
    lam = min(1., 2 * lam)
    while True:
      r = residuals(xk + lam * dx)
      if r is not None and np.linalg.norm(lin.step(*r)) <= (1 - lam / 4) * np.linalg.norm(dx):
        break
      lam /= 2
      if lam < 1e-6: raise RuntimeError('Gauss-Newton: no progress')
    xk = xk + lam * dx
    print('%s iter %2d: |F1|^2 = %.6e, |F2| = %.1e, |dx| = %.1e, damping %.3g'
      % (label, it, np.sum(f1**2), np.linalg.norm(f2), np.linalg.norm(dx), lam))
    if lam == 1 and np.linalg.norm(dx) < tol * (1 + np.linalg.norm(xk)):
      return xk, it + 1, lin
  return xk, max_iter, lin

# Multiple shooting: nodes at t = 0, 1, ..., 5 and 5.5, state variables at the nodes;
# measured concentrations as initial guess for the states, all rate constants equal to one
nodes = [0, 2, 4, 6, 8, 10, 11]    # indices into t_meas
interval = ca.integrator('interval', 'rosenbrock', dae, 0, [0.5, 1.0], opts)
last = ca.integrator('last', 'rosenbrock', dae, 0, 0.5, opts)
n_nodes = len(nodes)


def ms_residuals(w):
  s = ca.reshape(w[:7 * n_nodes], 7, n_nodes)
  pp = w[7 * n_nodes:]
  pred = [s[:, 0]]
  match = []
  for k in range(n_nodes - 1):
    I = interval if nodes[k + 1] - nodes[k] == 2 else last
    xf = I(x0=s[:, k], p=pp)['xf']
    pred += ca.horzsplit(xf)
    match.append(xf[:, -1] - s[:, k + 1])
  return ca.vec(ca.horzcat(*pred) - data) / sigma, ca.vertcat(*match)

w0 = np.concatenate([data[:, nodes].ravel(order='F'), np.ones(11)])
w_ms, it_ms, lin = gauss_newton(lambda w: ms_residuals(w)[0], lambda w: ms_residuals(w)[1],
  w0, 'multiple shooting')
p_ms = w_ms[7 * n_nodes:]
# Standard errors from the linearized covariance, cf. (8.42)-(8.44)
std = np.sqrt(np.diag(lin.covariance()))[7 * n_nodes:]

# Single shooting: the initial state and the rate constants
full = ca.integrator('full', 'rosenbrock', dae, 0, t_meas[1:], opts)


def ss_residuals(w):
  pred = ca.horzcat(w[:7], full(x0=w[:7], p=w[7:])['xf'])
  return ca.vec(pred - data) / sigma

w0 = np.concatenate([data[:, 0], np.ones(11)])
try:
  w_ss, it_ss, _ = gauss_newton(ss_residuals, lambda w: ca.MX(0, 1), w0, 'single shooting')
except RuntimeError as e:
  w_ss, it_ss = None, None
  print('single shooting failed:', str(e).splitlines()[-1])

print('\nGauss-Newton iterations: multiple shooting %d, single shooting %s' % (it_ms, it_ss))
print('%4s %12s %12s %10s' % ('', 'true', 'estimate', 'std. error'))
for i in range(11):
  print('p%-3d %12.4g %12.4g %9.1f%%' % (i + 1, p_true[i], p_ms[i], 100 * std[i] / abs(p_ms[i])))
