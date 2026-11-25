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

public static class direct_collocation {
  // Polynomials as coefficient arrays, highest power first (like numpy.poly1d)
  static double[] polymul(double[] a, double[] b) {
    var r = new double[a.Length + b.Length - 1];
    for (int i = 0; i < a.Length; i++) for (int j = 0; j < b.Length; j++) r[i + j] += a[i] * b[j];
    return r;
  }

  static double polyval(double[] p, double x) {
    double y = 0;
    foreach (double c in p) y = y * x + c;
    return y;
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

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Degree of interpolating polynomial
    int d = 3;

    // Get collocation points
    double[] tau_root = ca.collocation_points(d, "legendre").Prepend(0.0).ToArray();

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
          p = polymul(p, [1, -tau_root[r]]);
          for (int i = 0; i < p.Length; i++) p[i] /= (tau_root[j] - tau_root[r]);
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

    // Formulate the NLP
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym($"U_{k}");
      w.Add(Uk);
      lbw.Add(-1);
      ubw.Add(1);
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
    MX w_all = ca.vertcat(w);
    MX g_all = ca.vertcat(g);
    MX x_plot_all = ca.horzcat(x_plot);
    MX u_plot_all = ca.horzcat(u_plot);

    // Create an NLP solver
    var prob = new Dictionary<string, MX> { ["f"] = J, ["x"] = w_all, ["g"] = g_all };
    Function solver = ca.nlpsol("solver", "ipopt", prob);

    // Function to get x and u trajectories from w
    var trajectories = ca.Function("trajectories", [w_all], [x_plot_all, u_plot_all],
        ["w"], ["x", "u"]);

    // Solve the NLP
    var sol = solver.call(new Dictionary<string, DM> {
      ["x0"] = w0.ToArray(), ["lbx"] = lbw.ToArray(), ["ubx"] = ubw.ToArray(),
      ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() });
    List<DM> traj = trajectories.call(new[] { sol["x"] });
    double[,] x_opt = traj[0].ToArray(); // to C# array
    double[,] u_opt = traj[1].ToArray(); // to C# array

    // Plot the result
    DM tgrid = ca.linspace(0.0, T, N + 1);
    // plot: x_opt[0, ..], x_opt[1, ..] and (step) u_opt[0, ..] versus tgrid
  }
}
