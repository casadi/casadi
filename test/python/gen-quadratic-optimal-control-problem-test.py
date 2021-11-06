from casadi import SX, vertcat, nlpsol
import numpy as np
from numpy import inf
from os.path import join
from os import makedirs
import argparse

parser = argparse.ArgumentParser(
    description=
    'Generate linear MPC test with quadratic cost without state constraints')
parser.add_argument('--nx',
                    required=True,
                    type=int,
                    help='Dimension of state vector x')
parser.add_argument('--nu',
                    required=True,
                    type=int,
                    help='Dimension of input vector u')
parser.add_argument('--seed',
                    required=True,
                    type=int,
                    help='Seed for the random number generator')
parser.add_argument('--output', required=True, type=str, help='Output folder')
parser.add_argument('--name', required=True, type=str, help='Name of the test')
args = parser.parse_args()

nx = args.nx
nu = args.nu
seed = args.seed
rng = np.random.default_rng(seed)

A = rng.normal(0, 1, (nx, nx))
B = rng.normal(0, 1, (nx, nu))

f = lambda x, u: A @ x + B @ u

u = SX.sym('u', nu, 1)
x = SX.sym('x', nx, 1)
unknowns = vertcat(u, x)

# Initial state and state after one time step
x0 = np.ones((nx, 1))

# Weight matrices Q (state) and R (input)
Q = 1e1 * np.eye(nx)
R = 1e0 * np.eye(nu)

# Objective
objective = u.T @ R @ u + x.T @ Q @ x

# Box constraints on input
u_lb, u_ub = -1, +1
bounds = {
    'lbx': [u_lb] * nu  # u
    + [-inf] * nx,  # x
    'ubx': [u_ub] * nu  # u
    + [inf] * nx,  # x
    'lbg': [0] * nx,
    'ubg': [0] * nx,
}

# Dynamics constraints
constraints = [
    f(x0, u) - x,
]

# Solve using IPOPT
cs_nlp = {
    'x': unknowns,
    'f': objective,
    'g': vertcat(*constraints),
}
cs_bounds = {
    'lbx': vertcat(*bounds['lbx']),
    'ubx': vertcat(*bounds['ubx']),
    'lbg': vertcat(*bounds['lbg']),
    'ubg': vertcat(*bounds['ubg']),
}
cs_opts = {
    'verbose': False,
    'ipopt.tol': 1e-5,
}

S = nlpsol('S', 'ipopt', cs_nlp, cs_opts)
r = S(x0=np.zeros((nu + nx, )), **cs_bounds)
lam = r['lam_g']
sol = r['x']
print(f'λ_g = {lam}')
print(f'sol: u = {sol[:nu]}')
print(f'     x = {sol[nu:]}')


def eigen(M, name=''):
    M = np.array(M)
    return name + ' << ' + ', '.join(map(str, M.reshape(-1))) + ';'


content = f"""#include <Eigen/Core>
#include <alpaqa/vec.hpp>

struct {args.name} {{

inline static constexpr unsigned nx() {{
    return {nx};
}}

inline static constexpr unsigned nu() {{
    return {nu};
}}

inline static alpaqa::mat A() {{
    auto A = alpaqa::mat({nx}, {nx});
    {eigen(A, 'A')}
    return A;
}}

inline static alpaqa::mat B() {{
    auto B = alpaqa::mat({nx}, {nu});
    {eigen(B, 'B')}
    return B;
}}

inline static alpaqa::vec solution() {{
    auto ux = alpaqa::vec({nx+nu});
    {eigen(sol, 'ux')}
    return ux;
}}

inline static alpaqa::vec lagrange_multipliers() {{
    auto lam = alpaqa::vec({nx});
    {eigen(lam, 'lam')}
    return lam;
}}

}};

#include <test-alm-gen.hpp>

TEST(ALM, {args.name}) {{
    do_test(make_problem<{args.name}>(), 
            {args.name}::solution(),
            {args.name}::lagrange_multipliers());
}}
"""

makedirs(args.output, exist_ok=True)
with open(join(args.output, args.name + '.cpp'), 'w') as f:
    f.write(content)
