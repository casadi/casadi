#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
"""Legendre-Gauss-Lobatto transcription of an optimal-control problem."""

import numpy as np
import casadi as cs
import matplotlib.pyplot as plt


def _legendre_all(N: int, x: np.ndarray):
    """
    Computes legendre polynomial P_N(x) and its derivative P'_N(x), and P''_N(x) using numerically
    stable three-term recurrence relations.
    """
    P0, P1 = np.ones_like(x), x.copy()
    dP0, dP1 = np.zeros_like(x), np.ones_like(x)
    d2P0, d2P1 = np.zeros_like(x), np.zeros_like(x)

    if N == 0:
        return P0, dP0, d2P0
    if N == 1:
        return P1, dP1, d2P1

    for k in range(1, N):
        P2 = ((2 * k + 1) * x * P1 - k * P0) / (k + 1)
        dP2 = ((2 * k + 1) * (P1 + x * dP1) - k * dP0) / (k + 1)
        d2P2 = ((2 * k + 1) * (2 * dP1 + x * d2P1) - k * d2P0) / (k + 1)

        P0, P1 = P1, P2
        dP0, dP1 = dP1, dP2
        d2P0, d2P1 = d2P1, d2P2

    return P1, dP1, d2P1


def lgl_nodes(N: int, tol: float = 1e-15, max_iter: int = 100) -> np.ndarray:
    """
    Given N intervals, computes N+1 Legendre-Gauss-Lobatto (LGL) points on [-1, 1]

    Roots of (1-x**2)P'_N(x) are LGL points

    Parameters:
        N (int): intervals/polynomial degree (returns N+1 grid points)
        tol (float): tolerance for newton ralphson method
        max_iter (int): maximum newton-raphson iterations.

    Returns:
        x (ndarray): a 1d array of N+1 LGL points from -1 to 1
    """
    if N < 1:
        raise ValueError("Degree N must be at least 1.")

    # inital guess
    k = np.arange(N + 1)
    x = -np.cos(np.pi * k / N)

    x[0] = -1.0
    x[-1] = 1.0

    # begin newton iterations
    if N > 1:
        x_int = x[1:-1].copy()
        for _ in range(max_iter):
            _, dP, d2P = _legendre_all(N, x_int)
            dx = dP / d2P
            x_int -= dx
            if np.max(np.abs(dx)) < tol:
                break
        x[1:-1] = x_int

    x[0] = -1.0
    x[-1] = 1.0

    # use symmetry to remove round-off errors
    x = 0.5 * (x - x[::-1])
    return x


def lgl_weights(N: int) -> np.ndarray:
    """
    Computes LGL quadrature weights for N+1 nodes on [-1, 1].

    wi=2/[(N)(N+1)(P_N(x_i))**2]

    Parameters:
        N (int): Polynomial degree (returns N+1 weights).

    Returns:
        w (ndarray): Array of shape (N+1,) containing integration weights.
    """
    if N < 1:
        raise ValueError("Degree N must be at least 1.")
    x = lgl_nodes(N)
    P_N, _, _ = _legendre_all(N, x)
    w = 2.0 / (N * (N + 1) * (P_N**2))
    # symmetry to avoid round-off errors
    w = 0.5 * (w + w[::-1])
    return w


def lgl_diff_matrix(N: int) -> np.ndarray:
    """
    Computes the (N+1)x(N+1) LGL differentiation matrix using
    exact polynomial ratios and the negative-sum trick.

    Derivative of the Lagrange interpolation polynomial passing through points
    (x1,y1),.....(x_(N+1),y_(N+1)), evaluated at x_1...x_(N+1) given by

    Yd=DY, where Y=[y1....y_(N+1)].T and Yd is the derivative of Y at LGL points.

    Parameters:
        N (int): Polynomial degree.

    Returns:
        D (ndarray): High-precision differentiation matrix of shape (N+1, N+1).
    """
    if N < 1:
        raise ValueError("Degree N must be at least 1.")

    x = lgl_nodes(N)
    P_N, _, _ = _legendre_all(N, x)

    P_ratio = np.outer(P_N, 1.0 / P_N)
    dX = x[:, None] - x[None, :]
    np.fill_diagonal(dX, 1.0)  # Avoid division by zero on diagonal
    D = P_ratio / dX
    np.fill_diagonal(D, 0.0)
    np.fill_diagonal(D, -np.sum(D, axis=1))
    D = 0.5 * (D - D[::-1, ::-1])
    return D


def lgl_setup(N: int):
    """
    Convenience wrapper returning LGL nodes (x), quadrature weights (w),
    and differentiation matrix (D).
    """
    x = lgl_nodes(N)
    w = lgl_weights(N)
    D = lgl_diff_matrix(N)
    return x, w, D


def solve_OCP(N: int, c=3):
    """
    template function to solve the OCP

    N is number of LGL intervals or degree of polynomial used to approximate the
    state

    c scaling factor to generate a denser grid

    returns dictionary of solution and data needed for plots
    """

    # OCP definition
    # Example 1 of "A Pseudospectral Method for the Optimal Control of
    # Constrained Feedback Linearizable Systems"
    # Qi Gong , Wei Kang, and I. Michael Ross
    # number of states, controls, dynamics, lagrange term, intial state
    # initial and final time
    nx = 2
    nu = 1
    x = cs.MX.sym("x", nx, 1)
    u = cs.MX.sym("u", nu, 1)
    xd = cs.Function("xd", [x, u], [cs.vertcat(x[1] ** 3, u)])
    lag = cs.Function("lag", [x, u], [4 * u**2])
    x_t0 = cs.DM([0, 1])
    t0 = 0
    tf = 2

    # decision variables
    # N+1 LGL points [-1,1] tau mapped to physical [t0,tf]
    # t=0.5(tf-t0)tau+0.5(tf+t0)
    # state and control at these points decision variables
    nlp = cs.Opti()
    X = nlp.variable(nx, N + 1)
    U = nlp.variable(nu, N + 1)
    x0 = X[:, 0]
    xf = X[:, -1]

    # obtain LGL nodes, quadrature weights, differentiation matrices
    tau, wi, D = lgl_setup(N)

    # define event, defect constraints
    event = cs.vertcat(x0 - x_t0)
    # differentiation matrix
    defect = (
        X @ D.T - (tf - t0) / 2 * xd.map(N + 1, "serial")(X, U) == 0
    )  # to extract dual vars later

    # intergral approx. with quadrature
    mayer = 4 * xf[0] + xf[1]
    lagrange = (
        0.5 * (tf - t0) * cs.dot(wi.reshape(1, -1), lag.map(N + 1, "serial")(X, U))
    )
    objective = lagrange + mayer

    nlp.minimize(objective)
    nlp.subject_to(defect)
    nlp.subject_to(event == 0)

    # ipopt solver and initial guess from linear interpolation
    nlp.solver("ipopt")
    nlp.set_initial(
        X, cs.vcat([cs.linspace(0, 1, N + 1).T, cs.linspace(0, 1, N + 1).T])
    )
    nlp.set_initial(U, cs.DM.ones(U.shape))
    sol = nlp.solve()

    # numerical state control state derivative time costates
    Xs = sol.value(X)
    Us = sol.value(U).reshape(nu, -1)
    Xds = sol.value(2 / (tf - t0) * X @ D.T)
    ts = (tf - t0) / 2 * tau + 0.5 * (tf + t0)

    # covector mapping principle (lagrange mult./quadrature weights)
    adjoint = sol.value(nlp.dual(defect) / cs.repmat(wi.reshape(1, -1), nx, 1))

    # hamiltonian for the OCP
    hamiltonian = lag.map(N + 1, "serial")(Xs, Us) + cs.sum1(adjoint * Xds)

    # Analytical reference for the local solution computed here.
    x1a = lambda t: -64 / (5 * (2 + t) ** 5) + 2 / 5
    x2a = lambda t: 4 / ((2 + t) ** 2)
    ua = lambda t: -8 / ((2 + t) ** 3)
    lam1a = lambda t: 4 * t**0
    lam2a = lambda t: 64 / ((2 + t) ** 3)
    td = (tf - t0) / 2 * lgl_nodes((N) * c) + 0.5 * (tf + t0)
    var_analytical = np.column_stack([x1a(td), x2a(td), ua(td), lam1a(td), lam2a(td)])
    var_numerical = np.hstack([Xs.T, Us.T, adjoint.T])
    label_a = [r"$x_{1a}$", r"$x_{2a}$", r"$u_a$", r"$\lambda_{1a}$", r"$\lambda_{2a}$"]
    label_n = [r"$x_1$", r"$x_2$", r"$u$", r"$\lambda_1$", r"$\lambda_2$"]

    data = {}
    data["var_analytical"] = var_analytical
    data["var_numerical"] = var_numerical
    data["label_a"] = label_a
    data["label_n"] = label_n
    data["td"] = td
    data["ts"] = ts
    data["adjoint"] = adjoint
    data["Xds"] = Xds
    data["Xs"] = Xs
    data["Us"] = Us
    data["hamiltonian"] = hamiltonian
    data["N"] = N
    data["constraint_jacobian_sparsity"] = cs.jacobian(nlp.g, nlp.x).sparsity()
    data["state_variable_count"] = X.numel()
    data["defect_constraint_count"] = defect.numel()

    return data


def plot_solution_verify(data):
    # primal dual solution vs analytical

    # numerical solution
    plt.figure()
    plt.plot(data["ts"], data["var_numerical"], "*", label=data["label_n"])
    plt.plot(data["td"], data["var_analytical"], "-", label=data["label_a"])
    plt.legend()
    plt.show()

    # Hamiltonian consistency check (a necessary optimality condition)
    # constant for time-invariant problem
    # zero for this specific problem
    plt.figure()
    plt.plot(
        data["ts"].flatten(),
        data["hamiltonian"].full().flatten(),
        "*",
        label="Numerical",
    )
    plt.plot(data["td"].flatten(), 0 * data["td"], "-", label="Analytical")
    plt.xlabel("time [s]")
    plt.ylabel("Hamiltonian")
    plt.legend()
    plt.axis("equal")
    plt.show()


def plot_jacobian_sparsity(data):
    """Show the dense state coupling from global polynomial differentiation."""
    sparsity = data["constraint_jacobian_sparsity"]
    rows, cols = sparsity.get_triplet()
    fig, ax = plt.subplots()
    ax.scatter(cols, rows, s=8, marker="s")
    # Opti stacks X before U and the dynamics before the initial conditions.
    ax.axvline(data["state_variable_count"] - 0.5, color="gray", linewidth=0.8)
    ax.axhline(data["defect_constraint_count"] - 0.5, color="gray", linewidth=0.8)
    ax.set_xlim(-0.5, sparsity.size2() - 0.5)
    ax.set_ylim(sparsity.size1() - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("Decision variable (states | controls)")
    ax.set_ylabel("Constraint (dynamics | initial state)")
    ax.set_title("Constraint Jacobian sparsity")
    fig.tight_layout()
    plt.show()


def spectral_convergence():
    """
    solve the given OCP for increasing grid sizes. A sudden dip in
    error with small grid sizes indicate spectral convergence
    """
    nmin = 5
    nmax = 50
    N = range(nmin, nmax, 5)
    var_err = []
    for i in N:
        data = solve_OCP(i, c=1)
        var_err.append(
            np.max(np.abs(data["var_analytical"] - data["var_numerical"]), axis=0)
        )
        print(var_err[-1].shape)
    print(type(var_err))
    plt.figure()
    plt.semilogy(np.array(N), np.row_stack(var_err), label=data["label_n"])
    plt.legend()
    plt.xlabel("Polynomial degree N (N+1 LGL points)")
    plt.ylabel("max. absolute error")
    plt.show()


if __name__ == "__main__":
    data = solve_OCP(N=25)
    plot_solution_verify(data)
    plot_jacobian_sparsity(data)
    spectral_convergence()
