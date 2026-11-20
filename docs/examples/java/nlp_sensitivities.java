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

public class nlp_sensitivities {
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

  static double[] concat(List<double[]> v) {
    return v.stream().flatMapToDouble(Arrays::stream).toArray();
  }

  static List<DM> zeros(Sparsity sp, int n) {
    List<DM> r = new ArrayList<>();
    for (int i = 0; i < n; i++) r.add(DM.zeros(sp));
    return r;
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("sqpmethod") || !ca.has_conic("qrqp")) {
      System.out.println("SKIP: nlpsol sqpmethod / conic qrqp not available");
      return;
    }
    double inf = Double.POSITIVE_INFINITY;

    // Degree of interpolating polynomial
    int d = 3;
    // Get collocation points
    double[] cp = ca.collocation_points(d, "legendre");
    double[] tau_root = new double[d + 1];
    System.arraycopy(cp, 0, tau_root, 1, d);
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
    SX xdot = ca.vertcat(List.of(ca.power(x2, 2).unaryMinus().plus(1).times(x1).minus(x2).plus(u), x1));

    // Objective term
    SX L = ca.power(x1, 2).plus(ca.power(x2, 2)).plus(ca.power(u, 2));

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

    // Perturb with P
    MX P = MX.sym("P", 2);
    Xk = Xk.plus(P);

    // Formulate the NLP
    for (int k = 0; k < N; k++) {

      // New NLP variable for the control
      MX Uk = MX.sym("U_" + k);
      w.add(Uk);
      lbw.add(new double[] {-1});
      ubw.add(new double[] {.85});
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
      MX Xk_end = Xk.times(D[0]);
      for (int j = 1; j < d + 1; j++) {
        // Expression for the state derivative at the collocation point
        MX xp = Xk.times(C[0][j]);
        for (int r = 0; r < d; r++) xp = xp.plus(Xc.get(r).times(C[r + 1][j]));

        // Append collocation equations
        @SuppressWarnings("unchecked")
        List<MX> fq = (List<MX>) f.invoke(Xc.get(j - 1), Uk);
        MX fj = fq.get(0), qj = fq.get(1);
        g.add(fj.times(h).minus(xp));
        lbg.add(new double[] {0, 0});
        ubg.add(new double[] {0, 0});

        // Add contribution to the end state
        Xk_end = Xk_end.plus(Xc.get(j - 1).times(D[j]));

        // Add contribution to quadrature function
        J = J.plus(qj.times(B[j] * h));
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
    MX wv = ca.vertcat(w);
    MX gv = ca.vertcat(g);
    MX x_plotv = ca.horzcat(x_plot);
    MX u_plotv = ca.horzcat(u_plot);
    DM w0v = new DM(concat(w0));
    DM lbwv = new DM(concat(lbw));
    DM ubwv = new DM(concat(ubw));
    DM lbgv = new DM(concat(lbg));
    DM ubgv = new DM(concat(ubg));

    // Create an NLP solver
    Map<String, MX> prob = Map.of("f", J, "x", wv, "g", gv, "p", P);

    // Function to get x and u trajectories from w
    Function trajectories = new Function("trajectories", List.of(wv), List.of(x_plotv, u_plotv),
        List.of("w"), List.of("x", "u"));

    // Create an NLP solver, using SQP and active-set QP for accurate multipliers
    Map<String, Object> opts = Map.of("qpsol", "qrqp",
        "qpsol_options", Map.of("print_iter", false, "error_on_fail", false), "print_time", false);
    Function solver = ca.nlpsol("solver", "sqpmethod", prob, opts);

    // Solve the NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", w0v, "lbx", lbwv, "ubx", ubwv,
        "lbg", lbgv, "ubg", ubgv, "p", 0));

    // Extract trajectories
    @SuppressWarnings("unchecked")
    List<DM> traj = (List<DM>) trajectories.invoke(sol.get("x"));
    double[][] x_opt = traj.get(0).full();
    double[][] u_opt = traj.get(1).full();

    // plot: x_opt[0], x_opt[1] and u_opt[0] (step) over tgrid
    DM tgrid = ca.linspace(new DM(0), new DM(T), N + 1);

    // High-level approach:
    // Use factory to e.g. to calculate Hessian of optimal f w.r.t. p
    Function hsolver = solver.factory("h", solver.name_in(), List.of("hess:f:p:p"));
    System.out.println("hsolver generated");
    Map<String, GenericMatrix> hsol = hsolver.invoke(Map.of("x0", sol.get("x"), "lam_x0", sol.get("lam_x"),
        "lam_g0", sol.get("lam_g"), "lbx", lbwv, "ubx", ubwv, "lbg", lbgv, "ubg", ubgv, "p", 0));
    System.out.println("Hessian of f w.r.t. p:");
    System.out.println(hsol.get("hess_f_p_p"));

    // Low-level approach: Calculate directional derivatives.
    // We will calculate two directions simultaneously
    int nfwd = 2;

    // Forward mode AD for the NLP solver object
    Function fwd_solver = solver.forward(nfwd);
    System.out.println("fwd_solver generated");

    // Seeds, initalized to zero
    List<DM> fwd_lbx = zeros(((DM) sol.get("x")).sparsity(), nfwd);
    List<DM> fwd_ubx = zeros(((DM) sol.get("x")).sparsity(), nfwd);
    List<DM> fwd_p = zeros(P.sparsity(), nfwd);
    List<DM> fwd_lbg = zeros(((DM) sol.get("g")).sparsity(), nfwd);
    List<DM> fwd_ubg = zeros(((DM) sol.get("g")).sparsity(), nfwd);

    // Let's preturb P
    fwd_p.get(0).set(0, 1); // first nonzero of P
    fwd_p.get(1).set(1, 1); // second nonzero of P

    // Calculate sensitivities using AD
    Map<String, Object> fwd_arg = new LinkedHashMap<>();
    fwd_arg.put("out_x", sol.get("x"));
    fwd_arg.put("out_lam_g", sol.get("lam_g"));
    fwd_arg.put("out_lam_x", sol.get("lam_x"));
    fwd_arg.put("out_f", sol.get("f"));
    fwd_arg.put("out_g", sol.get("g"));
    fwd_arg.put("lbx", lbwv);
    fwd_arg.put("ubx", ubwv);
    fwd_arg.put("lbg", lbgv);
    fwd_arg.put("ubg", ubgv);
    fwd_arg.put("fwd_lbx", ca.horzcat(fwd_lbx));
    fwd_arg.put("fwd_ubx", ca.horzcat(fwd_ubx));
    fwd_arg.put("fwd_lbg", ca.horzcat(fwd_lbg));
    fwd_arg.put("fwd_ubg", ca.horzcat(fwd_ubg));
    fwd_arg.put("p", 0);
    fwd_arg.put("fwd_p", ca.horzcat(fwd_p));
    Map<String, GenericMatrix> sol_fwd = fwd_solver.invoke(fwd_arg);

    // Calculate the same thing using finite differences
    h = 1e-3;
    List<Map<String, GenericMatrix>> pert = new ArrayList<>();
    for (int dd = 0; dd < nfwd; dd++) {
      DM dx = fwd_lbx.get(dd).plus(fwd_ubx.get(dd)).times(h);
      DM dg = fwd_lbg.get(dd).plus(fwd_ubg.get(dd)).times(h);
      pert.add(solver.invoke(Map.of("x0", sol.get("x"), "lam_g0", sol.get("lam_g"), "lam_x0", sol.get("lam_x"),
          "lbx", lbwv.plus(dx), "ubx", ubwv.plus(dx), "lbg", lbgv.plus(dg), "ubg", ubgv.plus(dg),
          "p", fwd_p.get(dd).times(h).plus(0))));
    }

    // Print the result
    DM solf = (DM) sol.get("f");
    for (String s : List.of("f")) {
      System.out.println("==========");
      System.out.println("Checking " + s);
      System.out.println("finite differences");
      for (int dd = 0; dd < nfwd; dd++) System.out.println(((DM) pert.get(dd).get(s)).minus(solf).div(h));
      System.out.println("AD fwd");
      DM M = (DM) sol_fwd.get("fwd_" + s);
      for (int dd = 0; dd < nfwd; dd++) System.out.println(Arrays.toString(M.get(new Slice(), dd).nonzeros()));
    }

    // Perturb again, in the opposite direction for second order derivatives
    List<Map<String, GenericMatrix>> pert2 = new ArrayList<>();
    for (int dd = 0; dd < nfwd; dd++) {
      DM dx = fwd_lbx.get(dd).plus(fwd_ubx.get(dd)).times(h);
      DM dg = fwd_lbg.get(dd).plus(fwd_ubg.get(dd)).times(h);
      pert2.add(solver.invoke(Map.of("x0", sol.get("x"), "lam_g0", sol.get("lam_g"), "lam_x0", sol.get("lam_x"),
          "lbx", lbwv.minus(dx), "ubx", ubwv.minus(dx), "lbg", lbgv.minus(dg), "ubg", ubgv.minus(dg),
          "p", fwd_p.get(dd).times(h).unaryMinus().plus(0))));
    }

    // Print the result
    for (String s : List.of("f")) {
      System.out.println("finite differences, second order: " + s);
      for (int dd = 0; dd < nfwd; dd++)
        System.out.println(((DM) pert.get(dd).get(s)).minus(solf.times(2)).plus((DM) pert2.get(dd).get(s)).div(h * h));
    }

    // Reverse mode AD for the NLP solver object
    int nadj = 1;
    Function adj_solver = solver.reverse(nadj);
    System.out.println("adj_solver generated");

    // Seeds, initalized to zero
    List<DM> adj_f = zeros(((DM) sol.get("f")).sparsity(), nadj);
    List<DM> adj_g = zeros(((DM) sol.get("g")).sparsity(), nadj);
    List<DM> adj_x = zeros(((DM) sol.get("x")).sparsity(), nadj);

    // Study which inputs influence f
    adj_f.get(0).set(0, 1);

    // Calculate sensitivities using AD
    Map<String, Object> adj_arg = new LinkedHashMap<>();
    adj_arg.put("out_x", sol.get("x"));
    adj_arg.put("out_lam_g", sol.get("lam_g"));
    adj_arg.put("out_lam_x", sol.get("lam_x"));
    adj_arg.put("out_f", sol.get("f"));
    adj_arg.put("out_g", sol.get("g"));
    adj_arg.put("lbx", lbwv);
    adj_arg.put("ubx", ubwv);
    adj_arg.put("lbg", lbgv);
    adj_arg.put("ubg", ubgv);
    adj_arg.put("adj_f", ca.horzcat(adj_f));
    adj_arg.put("adj_g", ca.horzcat(adj_g));
    adj_arg.put("p", 0);
    adj_arg.put("adj_x", ca.horzcat(adj_x));
    Map<String, GenericMatrix> sol_adj = adj_solver.invoke(adj_arg);

    // Print the result
    for (String s : List.of("p")) {
      System.out.println("==========");
      System.out.println("Checking " + s);
      System.out.println("Reverse mode AD");
      System.out.println(sol_adj.get("adj_" + s));
    }
  }
}
