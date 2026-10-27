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

public class vdp_collocation {
  // Polynomials as coefficient arrays, highest power first (numpy.poly1d convention)
  static double[] polymul(double[] a, double[] b) {
    double[] c = new double[a.length + b.length - 1];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < b.length; j++) c[i + j] += a[i] * b[j];
    return c;
  }

  static double polyval(double[] p, double x) {
    double r = 0;
    for (double c : p) r = r * x + c;
    return r;
  }

  static double[] polyder(double[] p) {
    int n = p.length - 1;
    double[] r = new double[Math.max(n, 1)];
    for (int i = 0; i < n; i++) r[i] = p[i] * (n - i);
    return r;
  }

  static double[] polyint(double[] p) {
    int n = p.length;
    double[] r = new double[n + 1];
    for (int i = 0; i < n; i++) r[i] = p[i] / (n - i);
    return r;
  }

  static boolean hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", Map.of("x", z, "f", z.times(z)),
        Map.of("ipopt.linear_solver", "ma27", "ipopt.print_level", 0, "print_time", false));
    s.invoke(Map.of("x0", 1));
    return (Boolean) s.stats().get("success");
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { System.out.println("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Degree of interpolating polynomial
    int d = 3;

    // Choose collocation points
    double[] cp = ca.collocation_points(d, "radau");
    double[] tau_root = new double[d + 1];
    System.arraycopy(cp, 0, tau_root, 1, d);

    // Coefficients of the collocation equation
    double[][] C = new double[d + 1][d + 1];

    // Coefficients of the continuity equation
    double[] D = new double[d + 1];

    // Coefficients of the quadrature function
    double[] F = new double[d + 1];

    // Construct polynomial basis
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      double[] p = {1};
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          double s = tau_root[j] - tau_root[r];
          p = polymul(p, new double[] {1 / s, -tau_root[r] / s});
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D[j] = polyval(p, 1.0);

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      double[] pder = polyder(p);
      for (int r = 0; r < d + 1; r++) {
        C[j][r] = polyval(pder, tau_root[r]);
      }

      // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      double[] pint = polyint(p);
      F[j] = polyval(pint, 1.0);
    }

    // Control discretization
    int nk = 20;

    // End time
    double tf = 10.0;

    // Size of the finite elements
    double h = tf / nk;

    // All collocation time points
    double[][] T = new double[nk][d + 1];
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < d + 1; j++)
        T[k][j] = h * (k + tau_root[j]);

    // Declare variables (use scalar graph)
    SX t = SX.sym("t");    // time
    SX u = SX.sym("u");    // control
    SX x = SX.sym("x", 2); // state

    // ODE rhs function and quadratures
    SX xdot = ca.vertcat(List.of(x.get(1).times(x.get(1)).unaryMinus().plus(1).times(x.get(0)).minus(x.get(1)).plus(u),
        x.get(0)));
    SX qdot = x.get(0).times(x.get(0)).plus(x.get(1).times(x.get(1))).plus(u.times(u));
    Function f = new Function("f", List.of(t, x, u), List.of(xdot, qdot));

    // Control bounds
    double u_min = -0.75;
    double u_max = 1.0;
    double u_init = 0.0;

    // State bounds and initial guess
    double[] x_min = {-inf, -inf};
    double[] x_max = { inf,  inf};
    double[] xi_min = {0.0, 1.0};
    double[] xi_max = {0.0, 1.0};
    double[] xf_min = {0.0, 0.0};
    double[] xf_max = {0.0, 0.0};
    double[] x_init = {0.0, 0.0};

    // Dimensions
    int nx = 2;
    int nu = 1;

    // Total number of variables
    int NX = nk * (d + 1) * nx; // Collocated states
    int NU = nk * nu;           // Parametrized controls
    int NXF = nx;               // Final state
    int NV = NX + NU + NXF;

    // NLP variable vector
    MX V = MX.sym("V", NV);

    // All variables with bounds and initial guess
    double[] vars_lb = new double[NV];
    double[] vars_ub = new double[NV];
    double[] vars_init = new double[NV];
    int offset = 0;

    // Get collocated states and parametrized control
    MX[][] X = new MX[nk + 1][d + 1];
    MX[] U = new MX[nk];
    for (int k = 0; k < nk; k++) {
      // Collocated states
      for (int j = 0; j < d + 1; j++) {
        // Get the expression for the state vector
        X[k][j] = V.get(new Slice(offset, offset + nx));

        // Add the initial condition
        System.arraycopy(x_init, 0, vars_init, offset, nx);

        // Add bounds
        if (k == 0 && j == 0) {
          System.arraycopy(xi_min, 0, vars_lb, offset, nx);
          System.arraycopy(xi_max, 0, vars_ub, offset, nx);
        } else {
          System.arraycopy(x_min, 0, vars_lb, offset, nx);
          System.arraycopy(x_max, 0, vars_ub, offset, nx);
        }
        offset += nx;
      }

      // Parametrized controls
      U[k] = V.get(new Slice(offset, offset + nu));
      Arrays.fill(vars_lb, offset, offset + nu, u_min);
      Arrays.fill(vars_ub, offset, offset + nu, u_max);
      Arrays.fill(vars_init, offset, offset + nu, u_init);
      offset += nu;
    }

    // State at end time
    X[nk][0] = V.get(new Slice(offset, offset + nx));
    System.arraycopy(xf_min, 0, vars_lb, offset, nx);
    System.arraycopy(xf_max, 0, vars_ub, offset, nx);
    System.arraycopy(x_init, 0, vars_init, offset, nx);
    offset += nx;

    // Constraint function for the NLP
    List<MX> g = new ArrayList<>();
    int ng = 0;

    // Objective function
    MX J = new MX(0);

    // For all finite elements
    for (int k = 0; k < nk; k++) {

      // For all collocation points
      for (int j = 1; j < d + 1; j++) {

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = new MX(0);
        for (int r = 0; r < d + 1; r++) {
          xp_jk = xp_jk.plus(X[k][r].times(C[r][j]));
        }

        // Add collocation equations to the NLP
        @SuppressWarnings("unchecked")
        List<MX> fq = (List<MX>) f.invoke(T[k][j], X[k][j], U[k]);
        MX fk = fq.get(0), qk = fq.get(1);
        g.add(fk.times(h).minus(xp_jk));
        ng += nx; // equality constraints

        // Add contribution to objective
        J = J.plus(qk.times(F[j] * h));
      }

      // Get an expression for the state at the end of the finite element
      MX xf_k = new MX(0);
      for (int r = 0; r < d + 1; r++) {
        xf_k = xf_k.plus(X[k][r].times(D[r]));
      }

      // Add continuity equation to NLP
      g.add(X[k + 1][0].minus(xf_k));
      ng += nx;
    }

    // Concatenate constraints
    MX gc = ca.vertcat(g);

    // NLP
    Map<String, MX> nlp = Map.of("x", V, "f", J, "g", gc);

    // ----
    // SOLVE THE NLP
    // ----

    // Set options
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    // opts.put("ipopt.max_iter", 4);
    opts.put("ipopt.linear_solver", "ma27");

    // Allocate an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, Object> arg = new LinkedHashMap<>();

    // Initial condition
    arg.put("x0", vars_init);

    // Bounds on x
    arg.put("lbx", vars_lb);
    arg.put("ubx", vars_ub);

    // Bounds on g
    arg.put("lbg", new double[ng]);
    arg.put("ubg", new double[ng]);

    // Solve the problem
    Map<String, GenericMatrix> res = solver.invoke(arg);

    // Print the optimal cost
    System.out.println("optimal cost:  " + ((DM) res.get("f")).to_double());

    // Retrieve the solution
    double[] v_opt = ((DM) res.get("x")).nonzeros();

    // Get values at the beginning of each finite element
    int stride = (d + 1) * nx + nu;
    double[] x0_opt = new double[nk + 1], x1_opt = new double[nk + 1], u_opt = new double[nk];
    for (int k = 0; k < nk + 1; k++) {
      x0_opt[k] = v_opt[k * stride];
      x1_opt[k] = v_opt[k * stride + 1];
    }
    for (int k = 0; k < nk; k++) u_opt[k] = v_opt[(d + 1) * nx + k * stride];
    DM tgrid = ca.linspace(new DM(0), new DM(tf), nk + 1);
    DM tgrid_u = ca.linspace(new DM(0), new DM(tf), nk);

    // plot: x0_opt, x1_opt over tgrid and u_opt (step) over tgrid_u, "Van der Pol optimization"
  }
}
