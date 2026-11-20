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

import org.casadi.*;
import java.util.*;

public class pseudospectral_collocation {

  /**
   * Computes legendre polynomial P_N(x) and its derivative P'_N(x), and P''_N(x) using numerically
   * stable three-term recurrence relations. Returns {P, dP, d2P}.
   */
  static double[][] _legendre_all(int N, double[] x) {
    int n = x.length;
    double[] P0 = new double[n], P1 = x.clone();
    double[] dP0 = new double[n], dP1 = new double[n];
    double[] d2P0 = new double[n], d2P1 = new double[n];
    Arrays.fill(P0, 1.0);
    Arrays.fill(dP1, 1.0);

    if (N == 0) return new double[][] {P0, dP0, d2P0};
    if (N == 1) return new double[][] {P1, dP1, d2P1};

    for (int k = 1; k < N; k++) {
      double[] P2 = new double[n], dP2 = new double[n], d2P2 = new double[n];
      for (int i = 0; i < n; i++) {
        P2[i] = ((2 * k + 1) * x[i] * P1[i] - k * P0[i]) / (k + 1);
        dP2[i] = ((2 * k + 1) * (P1[i] + x[i] * dP1[i]) - k * dP0[i]) / (k + 1);
        d2P2[i] = ((2 * k + 1) * (2 * dP1[i] + x[i] * d2P1[i]) - k * d2P0[i]) / (k + 1);
      }
      P0 = P1; P1 = P2;
      dP0 = dP1; dP1 = dP2;
      d2P0 = d2P1; d2P1 = d2P2;
    }

    return new double[][] {P1, dP1, d2P1};
  }

  static double[] lgl_nodes(int N) { return lgl_nodes(N, 1e-15, 100); }

  /**
   * Given N intervals, computes N+1 Legendre-Gauss-Lobatto (LGL) points on [-1, 1]
   *
   * Roots of (1-x**2)P'_N(x) are LGL points
   *
   * @param N intervals/polynomial degree (returns N+1 grid points)
   * @param tol tolerance for newton ralphson method
   * @param max_iter maximum newton-raphson iterations.
   * @return a 1d array of N+1 LGL points from -1 to 1
   */
  static double[] lgl_nodes(int N, double tol, int max_iter) {
    if (N < 1) throw new IllegalArgumentException("Degree N must be at least 1.");

    // inital guess
    double[] x = new double[N + 1];
    for (int k = 0; k <= N; k++) x[k] = -Math.cos(Math.PI * k / N);

    x[0] = -1.0;
    x[N] = 1.0;

    // begin newton iterations
    if (N > 1) {
      double[] x_int = Arrays.copyOfRange(x, 1, N);
      for (int it = 0; it < max_iter; it++) {
        double[][] l = _legendre_all(N, x_int);
        double maxdx = 0;
        for (int i = 0; i < x_int.length; i++) {
          double dx = l[1][i] / l[2][i];
          x_int[i] -= dx;
          maxdx = Math.max(maxdx, Math.abs(dx));
        }
        if (maxdx < tol) break;
      }
      System.arraycopy(x_int, 0, x, 1, N - 1);
    }

    x[0] = -1.0;
    x[N] = 1.0;

    // use symmetry to remove round-off errors
    double[] r = new double[N + 1];
    for (int i = 0; i <= N; i++) r[i] = 0.5 * (x[i] - x[N - i]);
    return r;
  }

  /**
   * Computes LGL quadrature weights for N+1 nodes on [-1, 1].
   *
   * wi=2/[(N)(N+1)(P_N(x_i))**2]
   */
  static double[] lgl_weights(int N) {
    if (N < 1) throw new IllegalArgumentException("Degree N must be at least 1.");
    double[] x = lgl_nodes(N);
    double[] P_N = _legendre_all(N, x)[0];
    double[] w = new double[N + 1];
    for (int i = 0; i <= N; i++) w[i] = 2.0 / (N * (N + 1) * (P_N[i] * P_N[i]));
    // symmetry to avoid round-off errors
    double[] r = new double[N + 1];
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
  static double[][] lgl_diff_matrix(int N) {
    if (N < 1) throw new IllegalArgumentException("Degree N must be at least 1.");

    double[] x = lgl_nodes(N);
    double[] P_N = _legendre_all(N, x)[0];

    int n = N + 1;
    double[][] D = new double[n][n];
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int j = 0; j < n; j++) {
        if (i == j) continue; // Avoid division by zero on diagonal
        D[i][j] = (P_N[i] / P_N[j]) / (x[i] - x[j]);
        s += D[i][j];
      }
      D[i][i] = -s;
    }
    double[][] r = new double[n][n];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) r[i][j] = 0.5 * (D[i][j] - D[n - 1 - i][n - 1 - j]);
    return r;
  }

  static class Data {
    double[][] var_analytical, var_numerical;
    List<String> label_a, label_n;
    double[] td, ts;
    DM adjoint, Xds, Xs, Us, hamiltonian;
    int N;
    Sparsity constraint_jacobian_sparsity;
    long state_variable_count, defect_constraint_count;
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
    Function xd = new Function("xd", List.of(x, u), List.of(ca.vertcat(List.of(ca.power(x.get(1), 3), u))));
    Function lag = new Function("lag", List.of(x, u), List.of(ca.power(u, 2).times(4)));
    DM x_t0 = new DM(new double[] {0, 1});
    double t0 = 0;
    double tf = 2;

    // decision variables
    // N+1 LGL points [-1,1] tau mapped to physical [t0,tf]
    // t=0.5(tf-t0)tau+0.5(tf+t0)
    // state and control at these points decision variables
    Opti nlp = new Opti();
    MX X = nlp.variable(nx, N + 1);
    MX U = nlp.variable(nu, N + 1);
    MX x0 = X.get(new Slice(), 0);
    MX xf = X.get(new Slice(), -1);

    // obtain LGL nodes, quadrature weights, differentiation matrices
    double[] tau = lgl_nodes(N);
    DM wi = new DM(new double[][] {lgl_weights(N)}); // row vector
    DM D = new DM(lgl_diff_matrix(N));

    // define event, defect constraints
    MX event = ca.vertcat(List.of(x0.minus(x_t0)));
    // differentiation matrix
    MX defect = ca.eq(ca.mtimes(X, D.T()).minus(((MX) xd.map(N + 1, "serial").invoke(X, U)).times((tf - t0) / 2)),
        0); // to extract dual vars later

    // intergral approx. with quadrature
    MX mayer = xf.get(0).times(4).plus(xf.get(1));
    MX lagrange = ca.dot(wi, (MX) lag.map(N + 1, "serial").invoke(X, U)).times(0.5 * (tf - t0));
    MX objective = lagrange.plus(mayer);

    nlp.minimize(objective);
    nlp.subject_to(defect);
    nlp.subject_to(ca.eq(event, 0));

    // ipopt solver and initial guess from linear interpolation
    nlp.solver("ipopt");
    DM lin = ca.linspace(new DM(0), new DM(1), N + 1).T();
    nlp.set_initial(X, ca.vertcat(List.of(lin, lin)));
    nlp.set_initial(U, DM.ones(U.size1(), U.size2()));
    OptiSol sol = nlp.solve();

    // numerical state control state derivative time costates
    DM Xs = sol.value(X);
    DM Us = ca.reshape(sol.value(U), nu, -1);
    DM Xds = sol.value(ca.mtimes(X.times(2 / (tf - t0)), D.T()));
    double[] ts = new double[N + 1];
    for (int i = 0; i <= N; i++) ts[i] = (tf - t0) / 2 * tau[i] + 0.5 * (tf + t0);

    // covector mapping principle (-Lagrange multipliers / quadrature weights)
    DM adjoint = sol.value(nlp.dual(defect).div(ca.repmat(wi, nx, 1))).unaryMinus();

    // hamiltonian for the OCP
    DM hamiltonian = ((DM) lag.map(N + 1, "serial").invoke(Xs, Us)).plus(ca.sum1(adjoint.times(Xds)));

    // Analytical reference for the local solution computed here.
    double[] tdn = lgl_nodes(N * c);
    double[] td = new double[tdn.length];
    double[][] var_analytical = new double[td.length][];
    for (int i = 0; i < td.length; i++) {
      double t = td[i] = (tf - t0) / 2 * tdn[i] + 0.5 * (tf + t0);
      double x1a = -64 / (5 * Math.pow(2 + t, 5)) + 2.0 / 5;
      double x2a = 4 / Math.pow(2 + t, 2);
      double ua = -8 / Math.pow(2 + t, 3);
      double lam1a = 4 * Math.pow(t, 0);
      double lam2a = 64 / Math.pow(2 + t, 3);
      var_analytical[i] = new double[] {x1a, x2a, ua, lam1a, lam2a};
    }
    double[][] var_numerical = ca.horzcat(List.of(Xs.T(), Us.T(), adjoint.T())).full();
    List<String> label_a = List.of("$x_{1a}$", "$x_{2a}$", "$u_a$", "$\\lambda_{1a}$", "$\\lambda_{2a}$");
    List<String> label_n = List.of("$x_1$", "$x_2$", "$u$", "$\\lambda_1$", "$\\lambda_2$");

    Data data = new Data();
    data.var_analytical = var_analytical;
    data.var_numerical = var_numerical;
    data.label_a = label_a;
    data.label_n = label_n;
    data.td = td;
    data.ts = ts;
    data.adjoint = adjoint;
    data.Xds = Xds;
    data.Xs = Xs;
    data.Us = Us;
    data.hamiltonian = hamiltonian;
    data.N = N;
    data.constraint_jacobian_sparsity = ca.jacobian(nlp.g(), nlp.x()).sparsity();
    data.state_variable_count = X.numel();
    data.defect_constraint_count = defect.numel();

    return data;
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
    List<long[]> triplet = sparsity.get_triplet(); // rows, cols
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
    List<double[]> var_err = new ArrayList<>();
    for (int i = nmin; i < nmax; i += 5) {
      Data data = solve_OCP(i, 1);
      int ncol = data.var_numerical[0].length;
      double[] err = new double[ncol];
      for (int r = 0; r < data.var_numerical.length; r++)
        for (int j = 0; j < ncol; j++)
          err[j] = Math.max(err[j], Math.abs(data.var_analytical[r][j] - data.var_numerical[r][j]));
      var_err.add(err);
      System.out.println("(" + var_err.get(var_err.size() - 1).length + ",)");
    }
    System.out.println(var_err.getClass());
    // plot: semilogy of var_err versus N, "max. absolute error"
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    Data data = solve_OCP(25, 3);
    plot_solution_verify(data);
    plot_jacobian_sparsity(data);
    spectral_convergence();
  }
}
