import numpy as np
import matplotlib.pyplot as plt
import yaml
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 15,
})

show_ψ = True

# xlim, ylim, zlim = [0, 5], [0, 5], [-1, 40]
xlim, ylim, zlim = [-5, 5], [-5, 5], [-1, 40]
levels = np.arange(0, 200, 5)

if show_ψ:
    levels = np.arange(0, 200, 5)



z_ub = np.array([-1, 2])
z_lb = np.array([-np.inf, -np.inf])

def meshplot(f, xlim, ylim, delta):
    x = np.arange(xlim[0], xlim[1], delta)
    y = np.arange(ylim[0], ylim[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
    return X, Y, Z


def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def g1(x, y):
    return -4 * x**2 + 0.25 * x**2 * y**2


def g2(x, y):
    return 0.125 * x**4 - y * x

def ψ(x, y, λ, Σ):
    tile = lambda v: np.tile(v.reshape(v.shape + (1, 1)), (1,) + x.shape)
    ζ = np.array([g1(x, y), g2(x, y)]) + tile(λ / Σ)
    ẑ = np.fmax(np.fmin(ζ, tile(z_ub)), tile(z_lb))
    d = ζ - ẑ
    ŷ = tile(Σ) * d
    return f(x, y) + 0.5 * (d[0] * ŷ[0] + d[1] * ŷ[1])

X, Y, Z = meshplot(f, xlim, ylim, 0.025)
Xg1, Yg1, Zg1 = meshplot(g1, xlim, ylim, 0.025)
Xg2, Yg2, Zg2 = meshplot(g2, xlim, ylim, 0.025)

with open("/tmp/1.yaml", "r") as file:
    data1 = yaml.safe_load(file)
with open("/tmp/2.yaml", "r") as file:
    data2 = yaml.safe_load(file)
with open("/tmp/3.yaml", "r") as file:
    data3 = yaml.safe_load(file)
data = [data1, data3, data2]
# data = [data1, data2]
alm_its = max(map(len, data))

fig, axs = plt.subplots(len(data), alm_its, sharex='all', sharey='all', figsize=(1.5*16,1.5*9))

for d in range(len(data)):
    for it in range(alm_its):
        if it >= len(data[d]): 
            axs[d][it].set_axis_off()
            continue
        x = np.array([e["x"] for e in data[d][it]["PANOC"]])
        λ = np.array(data[d][it]["y"])
        Σ = np.array(data[d][it]["Σ"])
        if show_ψ: 
            X, Y, Z = meshplot(lambda x, y: ψ(x, y, λ, Σ), xlim, ylim, 0.025)
        axs[d][it].contour(X, Y, Z, levels)
        axs[d][it].contour(Xg1, Yg1, Zg1, [z_ub[0]], colors=['black'])
        axs[d][it].contour(Xg2, Yg2, Zg2, [z_ub[1]], colors=['black'])
        axs[d][it].set_xlim(xlim)
        axs[d][it].set_ylim(ylim)
        axs[d][it].plot(x[:, 0], x[:, 1], 'r.-')
        # axs[d][it].set_title(f"inner it: {x.shape[0]}, $\Sigma=$[{Σ[0]:.2f}, {Σ[1]:.2f}]")
        axs[d][it].set_title(f"inner it: {x.shape[0]}")
        # if d == 0:
        #     x = np.array([e["x"] for e in data[2][it]["PANOC"]])
        #     axs[d][it].plot(x[:, 0], x[:, 1], 'g.--')
        if it == 0:
            axs[d][it].set_ylabel(['PANOC + L-BFGS', '2nd order PANOC L-BFGS', '2nd order PANOC Newton'][d])
        axs[d][it].set_aspect('equal')

plt.tight_layout()
plt.savefig("poster-alm-all.pdf")
plt.show()

selection = [0, 1, 2, 7]
nr, nc = 2, len(selection)//2
fig, axs = plt.subplots(nr, nc, sharex='all', sharey='all', figsize=(nc*4, nr*4))

d = 2
for i, it in enumerate(selection):
    a = axs[i//nc][i%nc]
    x = np.array([e["x"] for e in data[d][it]["PANOC"]])
    λ = np.array(data[d][it]["y"])
    Σ = np.array(data[d][it]["Σ"])
    X, Y, Z = meshplot(lambda x, y: ψ(x, y, λ, Σ), xlim, ylim, 0.025)
    # if it > 0:
    # X, Y, Z = meshplot(f, xlim, ylim, 0.025)
    a.contour(X, Y, Z, levels)
    # a.contour(Xg1, Yg1, Zg1, [z_ub[0]], colors=['black'])
    a.contour(Xg2, Yg2, Zg2, [z_ub[1]], colors=['black'])
    a.set_xlim(xlim)
    a.set_ylim(ylim)
    a.plot(x[:, 0], x[:, 1], 'r.-')
    a.set_title(f"Penalty {Σ[1]:.2f}")
    # axs[d][it].set_title(f"inner it: {x.shape[0]}, $\Sigma=$[{Σ[0]:.2f}, {Σ[1]:.2f}]")
    # axs[d][it].set_title(f"inner it: {x.shape[0]}")
    # if d == 0:
    #     x = np.array([e["x"] for e in data[2][it]["PANOC"]])
    #     axs[d][it].plot(x[:, 0], x[:, 1], 'g.--')
    # if it == 0:
    #     axs[d][it].set_ylabel(['PANOC + L-BFGS', '2nd order PANOC L-BFGS', '2nd order PANOC Newton'][d])
    a.set_aspect('equal')
plt.tight_layout()
plt.savefig("poster-alm.pdf")
plt.show()