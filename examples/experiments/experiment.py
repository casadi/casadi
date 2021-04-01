import yaml
import matplotlib.pyplot as plt
import numpy as np

with open("/tmp/panoc-experiment.yaml", "r") as f:
    data = yaml.safe_load(f)

it = 0
grad_ψ = np.array([e["∇ψ"] for e in data[it]["PANOC"]])
grad_ψ_hat = np.array([e["∇ψ_hat"] for e in data[it]["PANOC"]])
x = np.array([e["x"] for e in data[it]["PANOC"]])
x_hat = np.array([e["x_hat"] for e in data[it]["PANOC"]])
p = np.array([e["p"] for e in data[it]["PANOC"]])
L = np.array([e["L"] for e in data[it]["PANOC"]])
γ = np.array([e["γ"] for e in data[it]["PANOC"]])
ψ = np.array([e["ψ"] for e in data[it]["PANOC"]])
ε = np.array([e["ε"] for e in data[it]["PANOC"]])

from pprint import pprint

# pprint(data[it]["PANOC"])

fig, axs = plt.subplots(4, 2, sharex='all')

axs[0, 0].plot(ψ, '.-')
axs[0, 0].set_title("ψ")

axs[1, 0].plot(grad_ψ, '.-')
axs[1, 0].set_title("∇ψ")

axs[2, 0].plot(grad_ψ_hat, '.-')
axs[2, 0].set_title("∇ψ_hat")

axs[3, 0].plot(x - x[-1,:], '.-')
axs[3, 0].set_title("x")

axs[0, 1].semilogy(γ, '.-', label='γ')
axs[0, 1].semilogy(1. / L, '.-', label='1/L')
axs[0, 1].legend()
axs[0, 1].set_title("γ")

axs[1, 1].semilogy(ε, '.-')
axs[1, 1].set_title("ε")

axs[2, 1].plot(p, '.-')
axs[2, 1].set_title("p")

axs[3, 1].plot(x_hat - x_hat[-1,:], '.-')
axs[3, 1].set_title("x_hat")

fig, axs = plt.subplots(2, 2, sharex='all')
axs[0, 0].semilogy([e["Σ"] for e in data.values()], '.-')
axs[0, 0].set_title("Σ")

axs[1, 0].plot([e["y"] for e in data.values()], '.-')
axs[1, 0].set_title("y")

plt.show()