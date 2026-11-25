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

public static class nlp_sensitivities {
  // Polynomials as coefficient arrays, highest power first (numpy.poly1d convention)
  static double[] polymul(double[] a, double[] b) {
    var c = new double[a.Length + b.Length - 1];
    for (int i = 0; i < a.Length; i++)
      for (int j = 0; j < b.Length; j++) c[i + j] += a[i] * b[j];
    return c;
  }

  static double polyval(double[] p, double x) {
    double r = 0;
    foreach (double c in p) r = r * x + c;
    return r;
  }

  static double[] polyder(double[] p) {
    int n = p.Length - 1;
    var r = new double[Math.Max(n, 1)];
    for (int i = 0; i < n; i++) r[i] = p[i] * (n - i);
    return r;
  }

  static double[] polyint(double[] p) {
    int n = p.Length;
    var r = new double[n + 1];
    for (int i = 0; i < n; i++) r[i] = p[i] / (n - i);
    return r;
  }

  static List<DM> zeros(Sparsity sp, int n) {
    var r = new List<DM>();
    for (int i = 0; i < n; i++) r.Add(DM.zeros(sp));
    return r;
  }

  public static void Main() {
    if (!ca.has_nlpsol("sqpmethod") || !ca.has_conic("qrqp")) {
      Console.WriteLine("SKIP: nlpsol sqpmethod / conic qrqp not available");
      return;
    }

    // Degree of interpolating polynomial
    int d = 3;
    // Get collocation points
    double[] tau_root = ca.collocation_points(d, "legendre").Prepend(0).ToArray();
    // Coefficients of the collocation equation
    var C = new double[d + 1, d + 1];
    // Coefficients of the continuity equation
    var D = new double[d + 1];
    // Coefficients of the quadrature function
    var B = new double[d + 1];
    // Construct polynomial basis
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      double[] p = { 1 };
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          double s = tau_root[j] - tau_root[r];
          p = polymul(p, [1 / s, -tau_root[r] / s]);
        }
      }
      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D[j] = polyval(p, 1.0);
      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      double[] pder = polyder(p);
      for (int r = 0; r < d + 1; r++) {
        C[j, r] = polyval(pder, tau_root[r]);
      }

      // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      double[] pint = polyint(p);
      B[j] = polyval(pint, 1.0);
    }

    // Time horizon
    double T = 10.0;

    // Declare model variables
    SX x1 = SX.sym("x1");
    SX x2 = SX.sym("x2");
    SX x = ca.vertcat(x1, x2);
    SX u = SX.sym("u");

    // Model equations
    SX xdot = ca.vertcat((1 - ca.pow(x2, 2)) * x1 - x2 + u, x1);

    // Objective term
    SX L = ca.pow(x1, 2) + ca.pow(x2, 2) + ca.pow(u, 2);

    // Continuous time dynamics
    var f = ca.Function("f", [x, u], [xdot, L], ["x", "u"], ["xdot", "L"]);

    // Control discretization
    int N = 20; // number of control intervals
    double h = T / N;

    // Start with an empty NLP
    var w = new List<MX>();
    var w0 = new List<double>();
    var lbw = new List<double>();
    var ubw = new List<double>();
    MX J = 0;
    var g = new List<MX>();
    var lbg = new List<double>();
    var ubg = new List<double>();

    // For plotting x and u given w
    var x_plot = new List<MX>();
    var u_plot = new List<MX>();

    // "Lift" initial conditions
    MX Xk = MX.sym("X0", 2);
    w.Add(Xk);
    lbw.AddRange([0, 1]);
    ubw.AddRange([0, 1]);
    w0.AddRange([0, 1]);
    x_plot.Add(Xk);

    // Perturb with P
    MX P = MX.sym("P", 2);
    Xk += P;

    // Formulate the NLP
    for (int k = 0; k < N; k++) {

      // New NLP variable for the control
      MX Uk = MX.sym($"U_{k}");
      w.Add(Uk);
      lbw.Add(-1);
      ubw.Add(.85);
      w0.Add(0);
      u_plot.Add(Uk);

      // State at collocation points
      var Xc = new List<MX>();
      for (int j = 0; j < d; j++) {
        MX Xkj = MX.sym($"X_{k}_{j}", 2);
        Xc.Add(Xkj);
        w.Add(Xkj);
        lbw.AddRange([-0.25, -ca.inf]);
        ubw.AddRange([ca.inf, ca.inf]);
        w0.AddRange([0, 0]);
      }

      // Loop over collocation points
      MX Xk_end = D[0] * Xk;
      for (int j = 1; j < d + 1; j++) {
        // Expression for the state derivative at the collocation point
        MX xp = C[0, j] * Xk;
        for (int r = 0; r < d; r++) xp += C[r + 1, j] * Xc[r];

        // Append collocation equations
        List<MX> fq = f.call([Xc[j - 1], Uk]);
        MX fj = fq[0], qj = fq[1];
        g.Add(h * fj - xp);
        lbg.AddRange([0, 0]);
        ubg.AddRange([0, 0]);

        // Add contribution to the end state
        Xk_end += D[j] * Xc[j - 1];

        // Add contribution to quadrature function
        J += B[j] * qj * h;
      }

      // New NLP variable for state at end of interval
      Xk = MX.sym($"X_{k + 1}", 2);
      w.Add(Xk);
      lbw.AddRange([-0.25, -ca.inf]);
      ubw.AddRange([ca.inf, ca.inf]);
      w0.AddRange([0, 0]);
      x_plot.Add(Xk);

      // Add equality constraint
      g.Add(Xk_end - Xk);
      lbg.AddRange([0, 0]);
      ubg.AddRange([0, 0]);
    }

    // Concatenate vectors
    MX wv = ca.vertcat(w);
    MX gv = ca.vertcat(g);
    MX x_plotv = ca.horzcat(x_plot);
    MX u_plotv = ca.horzcat(u_plot);
    DM w0v = w0.ToArray();
    DM lbwv = lbw.ToArray();
    DM ubwv = ubw.ToArray();
    DM lbgv = lbg.ToArray();
    DM ubgv = ubg.ToArray();

    // Create an NLP solver
    var prob = new Dictionary<string, MX> { ["f"] = J, ["x"] = wv, ["g"] = gv, ["p"] = P };

    // Function to get x and u trajectories from w
    var trajectories = ca.Function("trajectories", [wv], [x_plotv, u_plotv],
        ["w"], ["x", "u"]);

    // Create an NLP solver, using SQP and active-set QP for accurate multipliers
    var opts = new Dictionary<string, object> {
      ["qpsol"] = "qrqp",
      ["qpsol_options"] = new Dictionary<string, object> { ["print_iter"] = false, ["error_on_fail"] = false },
      ["print_time"] = false };
    Function solver = ca.nlpsol("solver", "sqpmethod", prob, opts);

    // Solve the NLP
    var sol = solver.call(new Dictionary<string, DM> {
      ["x0"] = w0v, ["lbx"] = lbwv, ["ubx"] = ubwv,
      ["lbg"] = lbgv, ["ubg"] = ubgv, ["p"] = 0 });

    // Extract trajectories
    List<DM> traj = trajectories.call(new[] { sol["x"] });
    double[,] x_opt = traj[0].ToArray();
    double[,] u_opt = traj[1].ToArray();

    // plot: x_opt[0], x_opt[1] and u_opt[0] (step) over tgrid
    DM tgrid = ca.linspace(0, T, N + 1);

    // High-level approach:
    // Use factory to e.g. to calculate Hessian of optimal f w.r.t. p
    Function hsolver = solver.factory("h", solver.name_in(), ["hess:f:p:p"]);
    Console.WriteLine("hsolver generated");
    var hsol = hsolver.call(new Dictionary<string, DM> {
      ["x0"] = sol["x"], ["lam_x0"] = sol["lam_x"], ["lam_g0"] = sol["lam_g"],
      ["lbx"] = lbwv, ["ubx"] = ubwv, ["lbg"] = lbgv, ["ubg"] = ubgv, ["p"] = 0 });
    Console.WriteLine("Hessian of f w.r.t. p:");
    Console.WriteLine(hsol["hess_f_p_p"]);

    // Low-level approach: Calculate directional derivatives.
    // We will calculate two directions simultaneously
    int nfwd = 2;

    // Forward mode AD for the NLP solver object
    Function fwd_solver = solver.forward(nfwd);
    Console.WriteLine("fwd_solver generated");

    // Seeds, initalized to zero
    List<DM> fwd_lbx = zeros(sol["x"].sparsity(), nfwd);
    List<DM> fwd_ubx = zeros(sol["x"].sparsity(), nfwd);
    List<DM> fwd_p = zeros(P.sparsity(), nfwd);
    List<DM> fwd_lbg = zeros(sol["g"].sparsity(), nfwd);
    List<DM> fwd_ubg = zeros(sol["g"].sparsity(), nfwd);

    // Let's preturb P
    fwd_p[0][0] = 1; // first nonzero of P
    fwd_p[1][1] = 1; // second nonzero of P

    // Calculate sensitivities using AD
    var sol_fwd = fwd_solver.call(new Dictionary<string, DM> {
      ["out_x"] = sol["x"],
      ["out_lam_g"] = sol["lam_g"],
      ["out_lam_x"] = sol["lam_x"],
      ["out_f"] = sol["f"],
      ["out_g"] = sol["g"],
      ["lbx"] = lbwv,
      ["ubx"] = ubwv,
      ["lbg"] = lbgv,
      ["ubg"] = ubgv,
      ["fwd_lbx"] = ca.horzcat(fwd_lbx),
      ["fwd_ubx"] = ca.horzcat(fwd_ubx),
      ["fwd_lbg"] = ca.horzcat(fwd_lbg),
      ["fwd_ubg"] = ca.horzcat(fwd_ubg),
      ["p"] = 0,
      ["fwd_p"] = ca.horzcat(fwd_p) });

    // Calculate the same thing using finite differences
    h = 1e-3;
    var pert = new List<Dictionary<string, DM>>();
    for (int dd = 0; dd < nfwd; dd++) {
      DM dx = h * (fwd_lbx[dd] + fwd_ubx[dd]);
      DM dg = h * (fwd_lbg[dd] + fwd_ubg[dd]);
      pert.Add(solver.call(new Dictionary<string, DM> {
        ["x0"] = sol["x"], ["lam_g0"] = sol["lam_g"], ["lam_x0"] = sol["lam_x"],
        ["lbx"] = lbwv + dx, ["ubx"] = ubwv + dx, ["lbg"] = lbgv + dg, ["ubg"] = ubgv + dg,
        ["p"] = 0 + h * fwd_p[dd] }));
    }

    // Print the result
    foreach (string s in new[] { "f" }) {
      Console.WriteLine("==========");
      Console.WriteLine($"Checking {s}");
      Console.WriteLine("finite differences");
      for (int dd = 0; dd < nfwd; dd++) Console.WriteLine((pert[dd][s] - sol[s]) / h);
      Console.WriteLine("AD fwd");
      DM M = sol_fwd["fwd_" + s];
      for (int dd = 0; dd < nfwd; dd++) Console.WriteLine($"[{string.Join(", ", M[.., dd].nonzeros())}]");
    }

    // Perturb again, in the opposite direction for second order derivatives
    var pert2 = new List<Dictionary<string, DM>>();
    for (int dd = 0; dd < nfwd; dd++) {
      DM dx = h * (fwd_lbx[dd] + fwd_ubx[dd]);
      DM dg = h * (fwd_lbg[dd] + fwd_ubg[dd]);
      pert2.Add(solver.call(new Dictionary<string, DM> {
        ["x0"] = sol["x"], ["lam_g0"] = sol["lam_g"], ["lam_x0"] = sol["lam_x"],
        ["lbx"] = lbwv - dx, ["ubx"] = ubwv - dx, ["lbg"] = lbgv - dg, ["ubg"] = ubgv - dg,
        ["p"] = 0 - h * fwd_p[dd] }));
    }

    // Print the result
    foreach (string s in new[] { "f" }) {
      Console.WriteLine($"finite differences, second order: {s}");
      for (int dd = 0; dd < nfwd; dd++)
        Console.WriteLine((pert[dd][s] - 2 * sol[s] + pert2[dd][s]) / (h * h));
    }

    // Reverse mode AD for the NLP solver object
    int nadj = 1;
    Function adj_solver = solver.reverse(nadj);
    Console.WriteLine("adj_solver generated");

    // Seeds, initalized to zero
    List<DM> adj_f = zeros(sol["f"].sparsity(), nadj);
    List<DM> adj_g = zeros(sol["g"].sparsity(), nadj);
    List<DM> adj_x = zeros(sol["x"].sparsity(), nadj);

    // Study which inputs influence f
    adj_f[0][0] = 1;

    // Calculate sensitivities using AD
    var sol_adj = adj_solver.call(new Dictionary<string, DM> {
      ["out_x"] = sol["x"],
      ["out_lam_g"] = sol["lam_g"],
      ["out_lam_x"] = sol["lam_x"],
      ["out_f"] = sol["f"],
      ["out_g"] = sol["g"],
      ["lbx"] = lbwv,
      ["ubx"] = ubwv,
      ["lbg"] = lbgv,
      ["ubg"] = ubgv,
      ["adj_f"] = ca.horzcat(adj_f),
      ["adj_g"] = ca.horzcat(adj_g),
      ["p"] = 0,
      ["adj_x"] = ca.horzcat(adj_x) });

    // Print the result
    foreach (string s in new[] { "p" }) {
      Console.WriteLine("==========");
      Console.WriteLine($"Checking {s}");
      Console.WriteLine("Reverse mode AD");
      Console.WriteLine(sol_adj["adj_" + s]);
    }
  }
}
