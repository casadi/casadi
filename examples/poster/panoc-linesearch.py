import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects 
import yaml
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 15,
})

show_ψ = True

xlim, ylim = [2.95, 3.1], [1.875, 2.05]
# xlim, ylim = [2.85, 3.15], [1.95, 2.25]
levels = np.arange(-0.05, 0.51, 0.020)

boxy = 1.95

def meshplot(f, xlim, ylim, delta):
    x = np.arange(xlim[0], xlim[1], delta)
    y = np.arange(ylim[0], ylim[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
    return X, Y, Z


def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def grad_f(x, y):
    sq = lambda x: x**2
    g = np.zeros((2,))
    g[0] = 2 * (2 * x * (sq(x) + y - 11) + x + sq(y) - 7)
    g[1] = 2 * (sq(x) + 2 * y * (x + sq(y) - 7) + y - 11)
    return g

def hess_f(x, y):
    H = np.zeros((2,2))
    H[0,0] = 4 * (x**2 + y - 11) + 8 * x**2 + 2
    H[0,1] = 4 * x + 4 * y
    H[1,0] = 4 * x + 4 * y
    H[1,1] = 4 * (x + y**2 - 7) + 8 * y**2 + 2
    return H

X, Y, Z = meshplot(f, xlim, ylim, 0.0005)

def load_results(name, long = False):
    if long: name = "long-" + name
    with open("/tmp/" + name + ".yaml", "r") as file:
        return yaml.safe_load(file)


fig, ax = plt.subplots(1, 1, sharex='all', sharey='all', figsize=(5,5))
ax.contour(X, Y, Z, levels)
cg1 = ax.contour(X, Y, Y, [boxy], colors='black', linewidths=0.8 )
plt.setp(cg1.collections,
         path_effects=[patheffects.withTickedStroke(spacing=7, linewidth=0.8)])
ax.set_xlim(xlim)
ax.set_ylim(ylim)
x = np.array([3.09, 1.9])
γ = 0.011
x_gd = x - γ * grad_f(x[0], x[1])
x_nw = x - np.linalg.solve(hess_f(x[0], x[1]), grad_f(x[0], x[1]))

ax.plot([x_gd[0], x_nw[0]], [x_gd[1], x_nw[1]], '-', color='k', linewidth=0.5)
ax.plot([x[0], x_gd[0]], [x[1], x_gd[1]], '.-', color='r', label='Projected gradient', markersize=10)
ax.plot([x[0], x_nw[0]], [x[1], x_nw[1]], '.-', color='purple', label='Quasi-Newton', markersize=10)
ax.text(x_nw[0] - 0.015, x_nw[1] + 0.003, r'$\tau = 1$')
ax.text(x_gd[0] - 0.024, x_gd[1] - 0.0025, r'$\tau = 0$')
ax.text(x[0] - 0.004, x[1] - 0.01, r'$x^k$')
x_ls = (x_gd + x_nw) / 2
ax.text(x_ls[0] - 0.014, x_ls[1] - 0.000, r'$x^{k+1}$')
ax.plot(x_ls[0], x_ls[1], 'k.', markersize=10)
ax.set_aspect('equal')
ax.set_yticks([1.9, 1.95, 2, 2.05])
ax.set_xticks(np.linspace(*xlim, 4, endpoint=True))
# ax.hlines(boxy, *xlim, linestyles=':', colors='k', linewidth=1)

ax.plot(3, 2, 'kx')
ax.legend()
ax.set_title('Himmelblau with box constraint')
plt.tight_layout()
plt.savefig("panoc-linesearch.pdf")

fig, ax = plt.subplots(1, 1, sharex='all', sharey='all', figsize=(6.4, 5))

def fbe_(x, γ):
    gr = grad_f(x[0], x[1])
    xhat = x - γ * gr
    if (xhat[1] > boxy): xhat[1] = boxy
    p = xhat - x
    fbe = f(x[0], x[1]) + 1. / (2 * γ) * np.linalg.norm(p)**2 + np.dot(gr, p)
    return fbe
fbe = np.vectorize(lambda x: fbe_(x, γ), signature='(1)->()')
interp = np.vectorize(lambda τ: (1. - τ) * x_gd + τ * x_nw, signature='()->(1)')

τ = np.linspace(0, 1, 100)
xτ = interp(τ)
fbexτ = fbe(xτ)
fv = lambda x: f(x[0], x[1])
fxτ = np.vectorize(fv, signature='(1)->()')(xτ)
τboxy = (boxy - x_gd[1]) / (x_nw[1] - x_gd[1])
ax.plot(τ, fxτ, '--', color='tab:blue', linewidth=1, label=r'Objective $\psi(x)$')
ax.plot(list(τ[τ < τboxy]) + [τboxy], list(fxτ[τ < τboxy]) + [100], '-', color='tab:blue', linewidth=2, label=r'Objective and constraints $\psi(x) + \delta_C(x)$')
ax.plot(τ, fbexτ, color='darkorange', linewidth=2, label=r'Forward-backward envelope $\varphi_\gamma(x)$')
ax.set_xlim([0, 1])
ax.set_ylim([0, 0.15])
ax.set_yticks(np.linspace(0, 0.15, 4, endpoint=True))
ax.text(0+0.01, -0.02, r'$\leftarrow$ Gradient step')
ax.text(1-0.37, -0.02, r'Quasi-Newton step $\rightarrow$')
# ax.axvline(τboxy, 0, 1, linestyle=':', color='k')
plt.legend()
plt.title('Line search over the FBE')
plt.tight_layout()
ax.set_xlabel(r'$\tau$')
plt.savefig("panoc-fbe.pdf")

plt.show()