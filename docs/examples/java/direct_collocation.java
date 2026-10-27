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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;

public class direct_collocation {
  // Polynomials as coefficient arrays, highest power first (like numpy.poly1d)
  static double[] polymul(double[] a, double[] b) {
    double[] r = new double[a.length + b.length - 1];
    for (int i = 0; i < a.length; i++) for (int j = 0; j < b.length; j++) r[i + j] += a[i] * b[j];
    return r;
  }

  static double polyval(double[] p, double x) {
    double y = 0;
    for (double c : p) y = y * x + c;
    return y;
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

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Degree of interpolating polynomial
    int d = 3;

    // Get collocation points
    double[] tau_root = new double[d + 1];
    System.arraycopy(ca.collocation_points(d, "legendre"), 0, tau_root, 1, d);

    // Coefficients of the collocation equation
    double[][] C = new double[d + 1][d + 1];

    // Coefficients of the continuity equation
    double[] D = new double[d + 1];

    // Coefficients of the quadrature function
    double[] B = new double[d + 1];

    // Construct polynomial basis
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      double[] p = {1};
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          p = polymul(p, new double[] {1, -tau_root[r]});
          for (int i = 0; i < p.length; i++) p[i] /= (tau_root[j] - tau_root[r]);
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
      B[j] = polyval(pint, 1.0);
    }

    // Time horizon
    double T = 10.;

    // Declare model variables
    SX x1 = SX.sym("x1");
    SX x2 = SX.sym("x2");
    SX x = ca.vertcat(List.of(x1, x2));
    SX u = SX.sym("u");

    // Model equations
    SX xdot = ca.vertcat(List.of(ca.minus(1, x2.power(2)).times(x1).minus(x2).plus(u), x1));

    // Objective term
    SX L = x1.power(2).plus(x2.power(2)).plus(u.power(2));

    // Continuous time dynamics
    Function f = new Function("f", List.of(x, u), List.of(xdot, L), List.of("x", "u"), List.of("xdot", "L"));

    // Control discretization
    int N = 20; // number of control intervals
    double h = T / N;

    // Start with an empty NLP
    List<MX> w = new ArrayList<>();
    List<double[]> w0 = new ArrayList<>();
    List<double[]> lbw = new ArrayList<>();
    List<double[]> ubw = new ArrayList<>();
    MX J = new MX(0);
    List<MX> g = new ArrayList<>();
    List<double[]> lbg = new ArrayList<>();
    List<double[]> ubg = new ArrayList<>();

    // For plotting x and u given w
    List<MX> x_plot = new ArrayList<>();
    List<MX> u_plot = new ArrayList<>();

    // "Lift" initial conditions
    MX Xk = MX.sym("X0", 2);
    w.add(Xk);
    lbw.add(new double[] {0, 1});
    ubw.add(new double[] {0, 1});
    w0.add(new double[] {0, 1});
    x_plot.add(Xk);

    // Formulate the NLP
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym("U_" + k);
      w.add(Uk);
      lbw.add(new double[] {-1});
      ubw.add(new double[] {1});
      w0.add(new double[] {0});
      u_plot.add(Uk);

      // State at collocation points
      List<MX> Xc = new ArrayList<>();
      for (int j = 0; j < d; j++) {
        MX Xkj = MX.sym("X_" + k + "_" + j, 2);
        Xc.add(Xkj);
        w.add(Xkj);
        lbw.add(new double[] {-0.25, -inf});
        ubw.add(new double[] {inf, inf});
        w0.add(new double[] {0, 0});
      }

      // Loop over collocation points
      MX Xk_end = ca.times(D[0], Xk);
      for (int j = 1; j < d + 1; j++) {
        // Expression for the state derivative at the collocation point
        MX xp = ca.times(C[0][j], Xk);
        for (int r = 0; r < d; r++) xp = xp.plus(ca.times(C[r + 1][j], Xc.get(r)));

        // Append collocation equations
        List<MX> fq = f.call(List.of(Xc.get(j - 1), Uk));
        MX fj = fq.get(0), qj = fq.get(1);
        g.add(ca.times(h, fj).minus(xp));
        lbg.add(new double[] {0, 0});
        ubg.add(new double[] {0, 0});

        // Add contribution to the end state
        Xk_end = Xk_end.plus(ca.times(D[j], Xc.get(j - 1)));

        // Add contribution to quadrature function
        J = J.plus(ca.times(B[j], qj).times(h));
      }

      // New NLP variable for state at end of interval
      Xk = MX.sym("X_" + (k + 1), 2);
      w.add(Xk);
      lbw.add(new double[] {-0.25, -inf});
      ubw.add(new double[] {inf, inf});
      w0.add(new double[] {0, 0});
      x_plot.add(Xk);

      // Add equality constraint
      g.add(Xk_end.minus(Xk));
      lbg.add(new double[] {0, 0});
      ubg.add(new double[] {0, 0});
    }

    // Concatenate vectors
    MX w_all = ca.vertcat(w);
    MX g_all = ca.vertcat(g);
    MX x_plot_all = ca.horzcat(x_plot);
    MX u_plot_all = ca.horzcat(u_plot);

    // Create an NLP solver
    Map<String, MX> prob = Map.of("f", J, "x", w_all, "g", g_all);
    Function solver = ca.nlpsol("solver", "ipopt", prob);

    // Function to get x and u trajectories from w
    Function trajectories = new Function("trajectories", List.of(w_all), List.of(x_plot_all, u_plot_all),
        List.of("w"), List.of("x", "u"));

    // Solve the NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", concatenate(w0), "lbx", concatenate(lbw),
        "ubx", concatenate(ubw), "lbg", concatenate(lbg), "ubg", concatenate(ubg)));
    List<DM> traj = trajectories.call(List.of((DM) sol.get("x")));
    double[][] x_opt = traj.get(0).full(); // to Java array
    double[][] u_opt = traj.get(1).full(); // to Java array

    // Plot the result
    DM tgrid = ca.linspace(0, T, N + 1);
    // plot: x_opt[0], x_opt[1] and (step) u_opt[0] versus tgrid
  }

  // numpy.concatenate of a list of 1-D arrays
  static double[] concatenate(List<double[]> parts) {
    int n = 0;
    for (double[] a : parts) n += a.length;
    double[] r = new double[n];
    int k = 0;
    for (double[] a : parts) for (double v : a) r[k++] = v;
    return r;
  }
}
