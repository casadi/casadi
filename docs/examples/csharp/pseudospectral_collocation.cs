//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
using System;
using System.Collections.Generic;
using System.Linq;
using CasADi;

public static class pseudospectral_collocation {

  /**
   * Computes legendre polynomial P_N(x) and its derivative P'_N(x), and P''_N(x) using numerically
   * stable three-term recurrence relations. Returns (P, dP, d2P).
   */
  static (double[] P, double[] dP, double[] d2P) _legendre_all(int N, double[] x) {
    int n = x.Length;
    double[] P0 = new double[n], P1 = (double[])x.Clone();
    double[] dP0 = new double[n], dP1 = new double[n];
    double[] d2P0 = new double[n], d2P1 = new double[n];
    Array.Fill(P0, 1.0);
    Array.Fill(dP1, 1.0);

    if (N == 0) return (P0, dP0, d2P0);
    if (N == 1) return (P1, dP1, d2P1);

    for (int k = 1; k < N; k++) {
      double[] P2 = new double[n], dP2 = new double[n], d2P2 = new double[n];
      for (int i = 0; i < n; i++) {
        P2[i] = ((2 * k + 1) * x[i] * P1[i] - k * P0[i]) / (k + 1);
        dP2[i] = ((2 * k + 1) * (P1[i] + x[i] * dP1[i]) - k * dP0[i]) / (k + 1);
        d2P2[i] = ((2 * k + 1) * (2 * dP1[i] + x[i] * d2P1[i]) - k * d2P0[i]) / (k + 1);
      }
      (P0, P1) = (P1, P2);
      (dP0, dP1) = (dP1, dP2);
      (d2P0, d2P1) = (d2P1, d2P2);
    }

    return (P1, dP1, d2P1);
  }

  /**
   * Given N intervals, computes N+1 Legendre-Gauss-Lobatto (LGL) points on [-1, 1]
   *
   * Roots of (1-x**2)P'_N(x) are LGL points
   *
   * N: intervals/polynomial degree (returns N+1 grid points)
   * tol: tolerance for newton ralphson method
   * max_iter: maximum newton-raphson iterations.
   * returns a 1d array of N+1 LGL points from -1 to 1
   */
  static double[] lgl_nodes(int N, double tol = 1e-15, int max_iter = 100) {
    if (N < 1) throw new ArgumentException("Degree N must be at least 1.");

    // inital guess
    var x = new double[N + 1];
    for (int k = 0; k <= N; k++) x[k] = -Math.Cos(Math.PI * k / N);

    x[0] = -1.0;
    x[N] = 1.0;

    // begin newton iterations
    if (N > 1) {
      double[] x_int = x[1..N];
      for (int it = 0; it < max_iter; it++) {
        var l = _legendre_all(N, x_int);
        double maxdx = 0;
        for (int i = 0; i < x_int.Length; i++) {
          double dx = l.dP[i] / l.d2P[i];
          x_int[i] -= dx;
          maxdx = Math.Max(maxdx, Math.Abs(dx));
        }
        if (maxdx < tol) break;
      }
      Array.Copy(x_int, 0, x, 1, N - 1);
    }

    x[0] = -1.0;
    x[N] = 1.0;

    // use symmetry to remove round-off errors
    var r = new double[N + 1];
    for (int i = 0; i <= N; i++) r[i] = 0.5 * (x[i] - x[N - i]);
    return r;
  }

  /**
   * Computes LGL quadrature weights for N+1 nodes on [-1, 1].
   *
   * wi=2/[(N)(N+1)(P_N(x_i))**2]
   */
  static double[] lgl_weights(int N) {
    if (N < 1) throw new ArgumentException("Degree N must be at least 1.");
    double[] x = lgl_nodes(N);
    double[] P_N = _legendre_all(N, x).P;
    var w = new double[N + 1];
    for (int i = 0; i <= N; i++) w[i] = 2.0 / (N * (N + 1) * (P_N[i] * P_N[i]));
    // symmetry to avoid round-off errors
    var r = new double[N + 1];
    for (int i = 0; i <= N; i++) r[i] = 0.5 * (w[i] + w[N - i]);
    return r;
  }

  /**
   * Computes the (N+1)x(N+1) LGL differentiation matrix using
   * exact polynomial ratios and the negative-sum trick.
   *
   * Derivative of the Lagrange interpolation polynomial passing through points
   * (x1,y1),.....(x_(N+1),y_(N+1)), evaluated at x_1...x_(N+1) given by
   *
   * Yd=DY, where Y=[y1....y_(N+1)].T and Yd is the derivative of Y at LGL points.
   */
  static double[,] lgl_diff_matrix(int N) {
    if (N < 1) throw new ArgumentException("Degree N must be at least 1.");

    double[] x = lgl_nodes(N);
    double[] P_N = _legendre_all(N, x).P;

    int n = N + 1;
    var D = new double[n, n];
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int j = 0; j < n; j++) {
        if (i == j) continue; // Avoid division by zero on diagonal
        D[i, j] = (P_N[i] / P_N[j]) / (x[i] - x[j]);
        s += D[i, j];
      }
      D[i, i] = -s;
    }
    var r = new double[n, n];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) r[i, j] = 0.5 * (D[i, j] - D[n - 1 - i, n - 1 - j]);
    return r;
  }

  class Data {
    public double[][] var_analytical;
    public double[,] var_numerical;
    public List<string> label_a, label_n;
    public double[] td, ts;
    public DM adjoint, Xds, Xs, Us, hamiltonian;
    public int N;
    public Sparsity constraint_jacobian_sparsity;
    public long state_variable_count, defect_constraint_count;
  }

  /**
   * template function to solve the OCP
   *
   * N is number of LGL intervals or degree of polynomial used to approximate the
   * state
   *
   * c scaling factor to generate a denser grid
   *
   * returns solution and data needed for plots
   */
  static Data solve_OCP(int N, int c) {

    // OCP definition
    // Example 1 of "A Pseudospectral Method for the Optimal Control of
    // Constrained Feedback Linearizable Systems"
    // Qi Gong , Wei Kang, and I. Michael Ross
    // number of states, controls, dynamics, lagrange term, intial state
    // initial and final time
    int nx = 2;
    int nu = 1;
    MX x = MX.sym("x", nx, 1);
    MX u = MX.sym("u", nu, 1);
    var xd = ca.Function("xd", [x, u], [ca.vertcat(ca.pow(x[1], 3), u)]);
    var lag = ca.Function("lag", [x, u], [4 * ca.pow(u, 2)]);
    DM x_t0 = new double[] { 0, 1 };
    double t0 = 0;
    double tf = 2;

    // decision variables
    // N+1 LGL points [-1,1] tau mapped to physical [t0,tf]
    // t=0.5(tf-t0)tau+0.5(tf+t0)
    // state and control at these points decision variables
    var nlp = ca.Opti();
    MX X = nlp.variable(nx, N + 1);
    MX U = nlp.variable(nu, N + 1);
    MX x0 = X[.., 0];
    MX xf = X[.., ^1];

    // obtain LGL nodes, quadrature weights, differentiation matrices
    double[] tau = lgl_nodes(N);
    DM wi = ((DM)lgl_weights(N)).T(); // row vector
    DM D = lgl_diff_matrix(N);

    // define event, defect constraints
    MX @event = ca.vertcat(x0 - x_t0);
    // differentiation matrix
    MX defect = ca.eq(ca.mtimes(X, D.T()) - (tf - t0) / 2 * xd.map(N + 1, "serial").Call(X, U),
        0); // to extract dual vars later

    // intergral approx. with quadrature
    MX mayer = 4 * xf[0] + xf[1];
    MX lagrange = 0.5 * (tf - t0) * ca.dot(wi, lag.map(N + 1, "serial").Call(X, U));
    MX objective = lagrange + mayer;

    nlp.minimize(objective);
    nlp.subject_to(defect);
    nlp.subject_to(ca.eq(@event, 0));

    // ipopt solver and initial guess from linear interpolation
    nlp.solver("ipopt");
    DM lin = ca.linspace(0, 1, N + 1).T();
    nlp.set_initial(X, ca.vertcat(lin, lin));
    nlp.set_initial(U, DM.ones(U.size1(), U.size2()));
    OptiSol sol = nlp.solve();

    // numerical state control state derivative time costates
    DM Xs = sol.value(X);
    DM Us = ca.reshape(sol.value(U), nu, -1);
    DM Xds = sol.value(ca.mtimes(2 / (tf - t0) * X, D.T()));
    var ts = new double[N + 1];
    for (int i = 0; i <= N; i++) ts[i] = (tf - t0) / 2 * tau[i] + 0.5 * (tf + t0);

    // covector mapping principle (-Lagrange multipliers / quadrature weights)
    DM adjoint = -sol.value(nlp.dual(defect) / ca.repmat(wi, nx, 1));

    // hamiltonian for the OCP
    DM hamiltonian = lag.map(N + 1, "serial").Call(Xs, Us) + ca.sum1(adjoint * Xds);

    // Analytical reference for the local solution computed here.
    double[] tdn = lgl_nodes(N * c);
    var td = new double[tdn.Length];
    var var_analytical = new double[td.Length][];
    for (int i = 0; i < td.Length; i++) {
      double t = td[i] = (tf - t0) / 2 * tdn[i] + 0.5 * (tf + t0);
      double x1a = -64 / (5 * Math.Pow(2 + t, 5)) + 2.0 / 5;
      double x2a = 4 / Math.Pow(2 + t, 2);
      double ua = -8 / Math.Pow(2 + t, 3);
      double lam1a = 4 * Math.Pow(t, 0);
      double lam2a = 64 / Math.Pow(2 + t, 3);
      var_analytical[i] = [x1a, x2a, ua, lam1a, lam2a];
    }
    double[,] var_numerical = ca.horzcat(Xs.T(), Us.T(), adjoint.T()).ToArray();
    var label_a = new List<string> { "$x_{1a}$", "$x_{2a}$", "$u_a$", "$\\lambda_{1a}$", "$\\lambda_{2a}$" };
    var label_n = new List<string> { "$x_1$", "$x_2$", "$u$", "$\\lambda_1$", "$\\lambda_2$" };

    return new Data {
      var_analytical = var_analytical,
      var_numerical = var_numerical,
      label_a = label_a,
      label_n = label_n,
      td = td,
      ts = ts,
      adjoint = adjoint,
      Xds = Xds,
      Xs = Xs,
      Us = Us,
      hamiltonian = hamiltonian,
      N = N,
      constraint_jacobian_sparsity = ca.jacobian(nlp.g, nlp.x).sparsity(),
      state_variable_count = X.numel(),
      defect_constraint_count = defect.numel(),
    };
  }

  static void plot_solution_verify(Data data) {
    // primal dual solution vs analytical

    // plot: data.var_numerical over data.ts ("*") vs data.var_analytical over data.td ("-")

    // Hamiltonian consistency check (a necessary optimality condition)
    // constant for time-invariant problem
    // zero for this specific problem
    // plot: data.hamiltonian over data.ts vs 0 over data.td ("Hamiltonian")
  }

  /** Show the dense state coupling from global polynomial differentiation. */
  static void plot_jacobian_sparsity(Data data) {
    Sparsity sparsity = data.constraint_jacobian_sparsity;
    var (rows, cols) = sparsity.get_triplet();
    // Opti stacks X before U and the dynamics before the initial conditions.
    // plot: scatter of the (cols, rows) triplet of the constraint Jacobian sparsity, with separators at
    //       data.state_variable_count and data.defect_constraint_count
  }

  /**
   * solve the given OCP for increasing grid sizes. A sudden dip in
   * error with small grid sizes indicate spectral convergence
   */
  static void spectral_convergence() {
    int nmin = 5;
    int nmax = 50;
    var var_err = new List<double[]>();
    for (int i = nmin; i < nmax; i += 5) {
      Data data = solve_OCP(i, 1);
      int ncol = data.var_numerical.GetLength(1);
      var err = new double[ncol];
      for (int r = 0; r < data.var_numerical.GetLength(0); r++)
        for (int j = 0; j < ncol; j++)
          err[j] = Math.Max(err[j], Math.Abs(data.var_analytical[r][j] - data.var_numerical[r, j]));
      var_err.Add(err);
      Console.WriteLine($"({var_err[^1].Length},)");
    }
    Console.WriteLine(var_err.GetType());
    // plot: semilogy of var_err versus N, "max. absolute error"
  }

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    Data data = solve_OCP(25, 3);
    plot_solution_verify(data);
    plot_jacobian_sparsity(data);
    spectral_convergence();
  }
}
