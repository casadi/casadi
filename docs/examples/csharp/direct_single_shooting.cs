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

#pragma warning disable CS0162 // unreachable code: the CVODES branch is switched off, as in the original

public static class direct_single_shooting {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    double T = 10.0; // Time horizon
    int N = 20; // number of control intervals

    // Declare model variables
    MX x1 = MX.sym("x1");
    MX x2 = MX.sym("x2");
    MX x = ca.vertcat(x1, x2);
    MX u = MX.sym("u");

    // Model equations
    MX xdot = ca.vertcat((1 - ca.pow(x2, 2)) * x1 - x2 + u, x1);

    // Objective term
    MX L = ca.pow(x1, 2) + ca.pow(x2, 2) + ca.pow(u, 2);

    // Formulate discrete time dynamics
    Function F;
    if (false) {
      // CVODES from the SUNDIALS suite
      var dae = new Dictionary<string, MX> { ["x"] = x, ["p"] = u, ["ode"] = xdot, ["quad"] = L };
      F = ca.integrator("F", "cvodes", dae, 0, T / N);
    } else {
      // Fixed step Runge-Kutta 4 integrator
      int M = 4; // RK4 steps per interval
      double DT = T / N / M;
      var f = ca.Function("f", [x, u], [xdot, L]);
      MX X0 = MX.sym("X0", 2);
      MX U = MX.sym("U");
      MX X = X0;
      MX Q = 0;
      for (int j = 0; j < M; j++) {
        List<MX> k1 = f.call([X, U]);
        List<MX> k2 = f.call([X + DT / 2 * k1[0], U]);
        List<MX> k3 = f.call([X + DT / 2 * k2[0], U]);
        List<MX> k4 = f.call([X + DT * k3[0], U]);
        X += DT / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
        Q += DT / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
      }
      F = ca.Function("F", [X0, U], [X, Q], ["x0", "p"], ["xf", "qf"]);
    }

    // Evaluate at a test point
    var Fk = F.call(new Dictionary<string, DM> { ["x0"] = new double[] { 0.2, 0.3 }, ["p"] = 0.4 });
    Console.WriteLine(Fk["xf"]);
    Console.WriteLine(Fk["qf"]);

    // Start with an empty NLP
    var w = new List<MX>();
    var w0 = new List<double>();
    var lbw = new List<double>();
    var ubw = new List<double>();
    MX J = 0;
    var g = new List<MX>();
    var lbg = new List<double>();
    var ubg = new List<double>();

    // Formulate the NLP
    MX Xk = new double[] { 0, 1 };
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym($"U_{k}");
      w.Add(Uk);
      lbw.Add(-1.0);
      ubw.Add(1.0);
      w0.Add(0.0);

      // Integrate till the end of the interval
      var Fkk = F.call(new Dictionary<string, MX> { ["x0"] = Xk, ["p"] = Uk });
      Xk = Fkk["xf"];
      J += Fkk["qf"];

      // Add inequality constraint
      g.Add(Xk[0]);
      lbg.Add(-0.25);
      ubg.Add(ca.inf);
    }

    // Create an NLP solver
    var prob = new Dictionary<string, MX> { ["f"] = J, ["x"] = ca.vertcat(w), ["g"] = ca.vertcat(g) };
    Function solver = ca.nlpsol("solver", "ipopt", prob);

    // Solve the NLP
    var sol = solver.call(new Dictionary<string, DM> {
      ["x0"] = w0.ToArray(), ["lbx"] = lbw.ToArray(), ["ubx"] = ubw.ToArray(),
      ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() });
    DM w_opt = sol["x"];

    // Plot the solution
    DM u_opt = w_opt;
    var x_opt = new List<DM> { new double[] { 0, 1 } };
    for (int k = 0; k < N; k++) {
      Fk = F.call(new Dictionary<string, DM> { ["x0"] = x_opt[^1], ["p"] = u_opt[k] });
      x_opt.Add(Fk["xf"]);
    }
    DM x1_opt = ca.vertcat(x_opt.Select(r => r[0]));
    DM x2_opt = ca.vertcat(x_opt.Select(r => r[1]));

    double[] tgrid = Enumerable.Range(0, N + 1).Select(k => T / N * k).ToArray();
    // plot: x1_opt (dashed), x2_opt over tgrid, and step plot of vertcat(nan, u_opt)
  }
}
