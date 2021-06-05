#%%
import casadi as cs
import numpy as np
import time
import matplotlib.pyplot as plt
import tikzplotlib
from bicycle_model import generate_problem

Ts = 0.05  # Sampling time
N_hor = 18  # Horizon length
N_sim = 60

multipleshooting = False
R_obstacle = 2

f, nlp, bounds, n_states, n_inputs, first_input_idx = generate_problem(Ts, N_hor, R_obstacle, multipleshooting)

# %%

import panocpy as pa
from tempfile import TemporaryDirectory
import os
from os.path import join

name = "mpcproblem"
f_prob = cs.Function("f", [nlp["x"], nlp["p"]], [nlp["f"]])
g_prob = cs.Function("g", [nlp["x"], nlp["p"]], [nlp["g"]])
cgen, n, m, num_p = pa.generate_casadi_problem(name, f_prob, g_prob)

with TemporaryDirectory(prefix="") as tmpdir:
    cfile = cgen.generate(tmpdir)
    sofile = join(tmpdir, f"{name}.so")
    os.system(f"cc -fPIC -shared -O3 -march=native {cfile} -o {sofile}")
    print(sofile)
    prob = pa.load_casadi_problem_with_param(sofile, n, m)

prob.C.lowerbound = bounds["lbx"]
prob.C.upperbound = bounds["ubx"]
prob.D.lowerbound = bounds["lbg"]
prob.D.upperbound = bounds["ubg"]

#%% PANOC params

lbfgsmem = 10
tol = 1e-18
verbose = False
fig_out_folder = "/home/pieter/Documents/KUL/2020-2021/Master-Thesis-Pieter-Pas/LaTeX/Thesis/fig"

from datetime import timedelta
pgaparams = {
    "max_iter": 1000,
    "print_interval": 1000 if verbose else 0,
    "stop_crit": pa.PANOCStopCrit.ApproxKKT,
}
gaapgaparams = pgaparams | {
    "limitedqr_mem" : lbfgsmem,
}
panocparams = pgaparams | {
    # "lbfgs_stepsize": pa.LBFGSStepsize.BasedOnGradientStepSize,
}
panoc_params_0 = panocparams | {
    "update_lipschitz_in_linesearch": False,
}
panoc_params_1 = panocparams | {
    "update_lipschitz_in_linesearch": True,
}

solvers = [
    pa.PGASolver(
        pa.PGAParams(**pgaparams),
    ),
    # pa.GAAPGASolver(
    #     pa.GAAPGAParams(**gaapgaparams),
    # ),
    pa.PANOCSolver(
        pa.PANOCParams(**panoc_params_0),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
    pa.PANOCSolver(
        pa.PANOCParams(**panoc_params_1),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
    pa.StructuredPANOCLBFGSSolver(
        pa.StructuredPANOCLBFGSParams(**panoc_params_0),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
    pa.StructuredPANOCLBFGSSolver(
        pa.StructuredPANOCLBFGSParams(**panoc_params_1),
        pa.LBFGSParams(memory=lbfgsmem),
    ),
]
almparams = pa.ALMParams(
    max_iter=10,
    max_time=timedelta(seconds=10),
    print_interval=1 if verbose else 0,
    preconditioning=False,
    ε=tol,
    δ=tol,
    ε_0=tol,
    Δ=10,
    Σ_0=(4e5 if multipleshooting else 4e5),
    max_total_num_retries=0,
    σ_0=1e-1,
    Σ_max=1e12,
    # single_penalty_factor=True,
)

#%%

resfpr = []
respsi = []
respsihat = []
resfbe = []
reslabels = ["PGA"]
# reslabels += ["GAAPGA"]
reslabels += ["PANOC (orig.)", "PANOC (impr.)"]
reslabels += ["2\\textsuperscript{nd} order PANOC (orig.)", "2\\textsuperscript{nd} order PANOC (impr.)"]

for reslabel, solver in zip(reslabels, solvers):
    almsolver = pa.ALMSolver(almparams, solver)
    y0 = np.zeros((m,))

    state = np.array([-5, 0, 0, 0])
    dest = np.array([5, 0.1, 0, 0])
    if multipleshooting:
        x_sol = np.concatenate((np.tile(state, N_hor), np.zeros((n_inputs * N_hor,))))
    else:
        x_sol = np.zeros((n,))
    assert x_sol.size == n
    y_sol = np.zeros((m,))


    def solve_ocp(state, y_sol, x_sol):
        state = np.reshape(state, (n_states,))
        prob.param = np.concatenate((state, dest))
        t0 = time.perf_counter()
        y_sol, x_sol, stats = almsolver(prob, y_sol, x_sol)
        t1 = time.perf_counter()
        return t1 - t0, stats, state, y_sol, x_sol

    class fpr_logger:
        def __init__(self):
            self.alm_it = -1
            self.fpr = []
            self.psi = []
            self.psi_hat = []
            self.x = []
            self.fbe = []
        
        def update(self, s):
            if s.k == 0:
                self.alm_it += 1
                self.fpr.append([])
                self.psi.append([])
                self.psi_hat.append([])
                self.x.append([])
                self.fbe.append([])
            self.fpr[self.alm_it].append(s.fpr)
            self.psi[self.alm_it].append(s.ψ)
            self.psi_hat[self.alm_it].append(s.ψ_hat)
            self.x[self.alm_it].append(np.linalg.norm(s.x))
            # try:
            #     self.fbe[self.alm_it].append(s.φγ)
            # except:
            φγ = 1 / (2 * s.γ) * s.norm_sq_p + np.dot(s.p, s.grad_ψ)
            self.fbe[self.alm_it].append(φγ)

    logger = fpr_logger()


    almsolver.inner_solver().set_progress_callback(logger.update)
    tdelta, stats, state, y_sol, x_sol = solve_ocp(state, y_sol, x_sol)
    print(reslabel)
    print(tdelta)
    print(stats['status'])
    resfpr.append(np.array(logger.fpr[0]))
    respsi.append(np.array(logger.psi[0]))
    respsihat.append(np.array(logger.psi_hat[0]))
    resfbe.append(np.array(logger.fbe[0]))

plt.figure()
for fpr, lbl in zip(resfpr, reslabels):
    # plt.semilogy(fpr, '-', label=f'{lbl} ({len(fpr)})')
    plt.semilogy(fpr, '-', label=f'{lbl})')
# plt.xlim([0, 200])
# plt.ylim([1e-4, 1e12])
plt.legend(loc='upper left')
plt.savefig(f'fpr-upd-ls-cond.pdf')
plt.title('FPR $\\|R_\\gamma(x)$\\|')
plt.xlabel('Iteration')
plt.ylabel('FPR')

plt.figure()
minpsi = min(map(np.min, respsihat))
print(f'{minpsi=}')
for psi, lbl in zip(respsi, reslabels):
    # plt.semilogy(psi - minpsi, '-', label=f'{lbl} ({len(psi)})')
    plt.semilogy(psi - minpsi, '-', label=f'{lbl}')
plt.xlim([0, pgaparams["max_iter"]])

plt.ylim([1e-10, 1e14])
plt.legend(loc='upper left', fontsize=13)
plt.savefig(f'psi-upd-ls-cond.pdf')
plt.title('Convergence of PGA and PANOC (Objective)')
plt.xlabel('Iteration')
plt.ylabel('$\\psi(x^k) - \\min_k \\psi(\\hat x^k)$')

fname = 'pga-panoc-ls-cond.tex'
tikzplotlib.save(join(fig_out_folder, fname), 
                    axis_width='0.8\linewidth',
                    extra_axis_parameters=['font=\small'])

plt.figure()
from pprint import pprint
# pprint(resfbe)
# minpsi = min(map(np.min, resfbe))
# print(f'{minpsi=}')
for psi, fbe, lbl in zip(respsi, resfbe, reslabels):
    # plt.semilogy(psi - minpsi, '-', label=f'{lbl} ({len(psi)})')
    plt.semilogy((psi - minpsi) + fbe, '-', label=f'{lbl}')
plt.xlim([0, pgaparams["max_iter"]])
plt.ylim([1e-10, 1e16])
# plt.legend(loc='upper left', fontsize=13)
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
plt.savefig(f'fbe-upd-ls-cond.pdf')
plt.title('Convergence of PGA and PANOC (FBE)')
plt.xlabel('Iteration')
plt.ylabel('$\\varphi_\\gamma(x^k) - \\min_k \\psi(\\hat x^k)$')

fname = 'pga-panoc-ls-cond-fbe.tex'
tikzplotlib.save(join(fig_out_folder, fname), 
                    axis_width='0.58\linewidth',
                    extra_axis_parameters=['font=\small'])

plt.show()

# plt.figure()
# for i, x in enumerate(logger.x):
#     plt.semilogy(np.array(x), '.-', label=f'ALM iter {i+1}')
# plt.xlim([0, 200])
# # plt.ylim([1e-11, 1e15])
# # plt.legend(loc='upper left')
# plt.savefig(f'x-upd-ls-cond={panocparams["update_lipschitz_in_linesearch"]}.pdf')
# plt.title('PANOC iterate norm $\\|x\\|$')
# plt.xlabel('Iteration')
# plt.ylabel('$\\|x\\|$')
# plt.show()

# %%
from pprint import pprint

greek = {
    'α': r'\(\alpha\)',
    'β': r'\(\beta\)',
    'γ': r'\(\gamma\)',
    'δ': r'\(\delta\)',
    'ϵ': r'\(\epsilon\)',
    'ε': r'\(\varepsilon\)',
    'ζ': r'\(\zeta\)',
    'η': r'\(\eta\)',
    'θ': r'\(\theta\)',
    'ι': r'\(\iota\)',
    'κ': r'\(\kappa\)',
    'λ': r'\(\lambda\)',
    'μ': r'\(\mu\)',
    'ν': r'\(\nu\)',
    'ξ': r'\(\xi\)',
    'ο': r'\(\omicron\)',
    'π': r'\(\pi\)',
    'ρ': r'\(\rho\)',
    'ς': r'\(\varsigma\)',
    'σ': r'\(\sigma\)',
    'τ': r'\(\tau\)',
    'υ': r'\(\upsilon\)',
    'φ': r'\(\phi\)',
    'χ': r'\(\chi\)',
    'ψ': r'\(\psi\)',
    'ω': r'\(\omega\)',
    'Α': r'\(\Alpha\)',
    'Β': r'\(\Beta\)',
    'Γ': r'\(\Gamma\)',
    'Δ': r'\(\Delta\)',
    'Ε': r'\(\Epsilon\)',
    'Ζ': r'\(\Zeta\)',
    'Η': r'\(\Eta\)',
    'Θ': r'\(\Theta\)',
    'Ι': r'\(\Iota\)',
    'Κ': r'\(\Kappa\)',
    'Λ': r'\(\Lambda\)',
    'Μ': r'\(\Mu\)',
    'Ν': r'\(\Nu\)',
    'Ξ': r'\(\Xi\)',
    'Ο': r'\(\Omicron\)',
    'Π': r'\(\Pi\)',
    'Ρ': r'\(\Rho\)',
    'Σ': r'\(\Varsigma\)',
    'Σ': r'\(\Sigma\)',
    'Τ': r'\(\Tau\)',
    'Υ': r'\(\Upsilon\)',
    'Φ': r'\(\Phi\)',
    'Χ': r'\(\Chi\)',
    'Ψ': r'\(\Psi\)',
    'Ω': r'\(\Omega\)',
    '_': r'\_',
}
greek_fix = lambda k: ''.join(map(lambda c: greek.get(c, c), k))

def print_settings_val(prefix, s):
    print(prefix + '} & ' + s + r' \\')


def print_settings(d, prefix = ''):
    if type(d) == dict:
        for k, v in d.items():
            kf = greek_fix(k)
            pf = prefix + '.' + kf if prefix else '\\texttt{' + kf
            print_settings(v, pf)
    elif type(d) == float:
        print_settings_val(prefix, f'{d:,.6g}')
    else:
        print_settings_val(prefix, str(d))

print_settings(        
    pa.PGAParams(**pgaparams).to_dict(),
)
print()
print_settings(        
    pa.PANOCParams(**panoc_params_0).to_dict(),
)
print()
print_settings(
    pa.StructuredPANOCLBFGSParams(**panoc_params_0).to_dict(),
)
print()
print_settings(        
    pa.LBFGSParams(memory=lbfgsmem).to_dict(),
)
print()
print_settings(        
    almparams.to_dict(),
)
# %%
