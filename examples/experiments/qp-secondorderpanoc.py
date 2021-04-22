import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects 
from LBFGS import LBFGS
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 15,
})


V = np.array([[1, -1], 
              [1, 1]])
A = 0.5 * V @ np.diag([1, 2]) @ V.T
print(V)
print(A)
print(np.linalg.eigvals(A))

x_lb = np.array([-np.inf, 0.25])
x_ub = np.array([0., +np.inf])
L = max(np.linalg.eigvals(A))
L = 8
γ = 0.95 * 1. / L

ψ = np.vectorize(lambda x: 0.5 * x.T @ A @ x, signature='(2)->()')
grad_ψ = np.vectorize(lambda x: A @ x, signature='(2)->(2)')
hess_ψ = np.vectorize(lambda x: A, signature='(2)->(2,2)')
Πc = np.vectorize(lambda x: np.minimum(np.maximum(x_lb, x), x_ub), signature='(2)->(2)')
Tγ = np.vectorize(lambda x: Πc(x - γ * grad_ψ(x)), signature='(2)->(2)')

def meshplot(f, xlim, ylim, delta):
    x = np.arange(xlim[0], xlim[1], delta)
    y = np.arange(ylim[0], ylim[1], delta)
    X, Y = np.meshgrid(x, y)
    XY = np.stack((X, Y), axis=2)
    Z = f(XY)
    return X, Y, Z

x0 = np.array([-2.4, 2])

xlim, ylim = [-2.5, 1], [-1, 2.5]
levels = np.linspace(0, 2.5**2*2, 50)

X, Y, Z = meshplot(ψ, xlim, ylim, 0.01)

fig, ax = plt.subplots(1, 1, figsize=(8,8))
ax.contour(X, Y, Z, levels)
if np.isfinite(x_lb[0]):
    cg0 = ax.contour(X, Y, -X, [-x_lb[0]], colors='black', linewidths=0.8, linestyles='-')
    plt.setp(cg0.collections,
            path_effects=[patheffects.withTickedStroke(spacing=7, linewidth=0.8)])
if np.isfinite(x_lb[1]):
    cg1 = ax.contour(X, Y, -Y, [-x_lb[1]], colors='black', linewidths=0.8, linestyles='-')
    plt.setp(cg1.collections,
            path_effects=[patheffects.withTickedStroke(spacing=7, linewidth=0.8)])
if np.isfinite(x_ub[0]):
    cg0 = ax.contour(X, Y, X, [x_ub[0]], colors='black', linewidths=0.8, linestyles='-')
    plt.setp(cg0.collections,
            path_effects=[patheffects.withTickedStroke(spacing=7, linewidth=0.8)])
if np.isfinite(x_ub[1]):
    cg1 = ax.contour(X, Y, Y, [x_ub[1]], colors='black', linewidths=0.8, linestyles='-')
    plt.setp(cg1.collections,
            path_effects=[patheffects.withTickedStroke(spacing=7, linewidth=0.8)])

N = 20
x_gd = np.nan * np.ones((N,2))
x_gd[0] = x0
for k in range(N-1):
    x_gd[k+1] = Tγ(x_gd[k])
plt.plot(x_gd[:,0], x_gd[:,1], 'r.-', label='PGA')

x_2nd = np.nan * np.ones((N,2))
x_2nd[0] = x0
τ = np.nan * np.ones((N,))
τ_min = 1. / 256
lbfgs = LBFGS(2, 5)
ε = 1e-12

for k in range(N-1):
    x = x_2nd[k]
    xhat = Tγ(x)
    p = xhat - x
    grad = grad_ψ(x)

    if np.linalg.norm((1 / γ) * p + (grad - grad_ψ(xhat)), ord=np.inf) < ε:
        x_2nd[k+1] = xhat
        x_2nd = x_2nd[:k+2]
        τ = τ[:k+1]
        break

    q = np.nan * np.ones((2,))
    J = []
    if k > 0:
        for i in range(2):
            gd = x[i] - γ * grad[i]
            if gd < x_lb[i] or x_ub[i] < gd:
                q[i] = p[i]
            else:
                J.append(i)
                q[i] = 0
        print(f'{J=}')
        Hq = hess_ψ(x) @ q
        q[J] = -grad[J] - Hq[J]
        if len(J) > 0:
            if lbfgs.apply(q, J) == False:
                q[J] *= γ
    τ[k] = 0. if k == 0 else 1.
    if len(J) == 0: τ[k] = 0

    # Line search
    σₖppγ = (1 - γ * L) * np.dot(p, p) / (2 * γ)
    φ = ψ(x) + 1 / (2 * γ) * np.dot(p, p) + np.dot(grad, p)
    while True:
        xkn = x + (1-τ[k]) * p + τ[k] * q
        xkn_hat = Tγ(xkn)
        pkn = xkn_hat - xkn 
        gradkn = grad_ψ(xkn)
        φkn = ψ(xkn) + 1 / (2 * γ) * np.dot(pkn, pkn) + np.dot(gradkn, pkn)
        ls_cond = φkn - (φ - σₖppγ)
        if ls_cond <= 0:
            break
        if τ[k] < τ_min:
            τ[k] = 0
            xkn = xhat
            xkn_hat = Tγ(xkn)
            pkn = xkn_hat - xkn 
            gradkn = grad_ψ(xkn)
            break
        τ[k] /= 2
    lbfgs.force_update(x, xkn, grad, gradkn)
    x_2nd[k+1] = xkn

plt.plot(x_2nd[:,0], x_2nd[:,1], 'b.-', linewidth=2, label=r'2\textsuperscript{nd} order PANOC + LBFGS (Python)')

import yaml
def load_results(name, long = False):
    if long: name = "long-" + name
    with open("/tmp/" + name + ".yaml", "r") as file:
        return yaml.safe_load(file)
cpp_data = load_results("panoc-2lbfgs-2")
x_cpp = np.array([e["x"] for e in cpp_data[0]["PANOC"]])
plt.plot(x_cpp[:,0], x_cpp[:,1], 'gx:', linewidth=2, markersize=8, label=r'2\textsuperscript{nd} order PANOC + LBFGS (C++)')
plt.legend()

plt.figure()
plt.plot(τ, '.-')
plt.ylabel(r'$\tau$')
plt.xlabel('Iteration')

plt.show()
