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

// Lotka-Volerra problem from mintoc.de
public static class lotka_volterra_minlp {
  public static void Main() {
    if (!ca.has_nlpsol("bonmin")) { Console.WriteLine("SKIP: nlpsol bonmin not available"); return; }

    double T = 12.0; // Time horizon
    int N = 40; // number of control intervals

    // Model parameters
    double c0 = 0.4;
    double c1 = 0.2;

    // Initial condition for x
    double[] x0 = { 0.5, 0.7 };

    // Bounds on x
    double[] lbx = { 0.0, 0.0 };
    double[] ubx = { 2.0, 2.0 };

    // Bounds on u
    double lbu = 0.0;
    double ubu = 1.0;

    // Declare model variables
    SX x = SX.sym("x", 2);
    SX u = SX.sym("u");

    // Model equations
    SX xdot = ca.vertcat(x[0] - x[0] * x[1] - c0 * x[0] * u,
                      -x[1] + x[0] * x[1] - c1 * x[1] * u);

    // Objective term
    SX L = ca.pow(x[0] - 1, 2) + ca.pow(x[1] - 1, 2) + 1e-4 * ca.pow(u, 2);

    // Formulate discrete time dynamics
    Function F;
    if (false) {
      // CVODES from the SUNDIALS suite
      var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = u, ["ode"] = xdot, ["quad"] = L };
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

    // Initial guess for u
    var u_start = Enumerable.Repeat((DM)0.0, N).ToList();

    // Get a feasible trajectory as an initial guess
    DM xk = x0;
    var x_start = new List<DM> { xk };
    for (int k = 0; k < N; k++) {
      xk = F.call(new Dictionary<string, DM> { ["x0"] = xk, ["p"] = u_start[k] })["xf"];
      x_start.Add(xk);
    }

    // Start with an empty NLP
    var w = new List<MX>();
    var w0 = new List<DM>();
    var lbw = new List<double>();
    var ubw = new List<double>();
    var discrete = new List<bool>();
    MX J = 0;
    var g = new List<MX>();
    var lbg = new List<double>();
    var ubg = new List<double>();

    // "Lift" initial conditions
    {
      MX X0 = MX.sym("X0", 2);
      w.Add(X0);
      lbw.AddRange(x0);
      ubw.AddRange(x0);
      w0.Add(x_start[0]);
      discrete.AddRange([false, false]);

      // Formulate the NLP
      MX Xk = X0;
      for (int k = 0; k < N; k++) {
        // New NLP variable for the control
        MX Uk = MX.sym($"U_{k}");
        w.Add(Uk);
        lbw.Add(lbu);
        ubw.Add(ubu);
        w0.Add(u_start[k]);
        discrete.Add(true);

        // Integrate till the end of the interval
        var Fk = F.call(new Dictionary<string, MX> { ["x0"] = Xk, ["p"] = Uk });
        MX Xk_end = Fk["xf"];
        J += Fk["qf"];

        // New NLP variable for state at end of interval
        Xk = MX.sym($"X_{k + 1}", 2);
        w.Add(Xk);
        lbw.AddRange(lbx);
        ubw.AddRange(ubx);
        w0.Add(x_start[k + 1]);
        discrete.AddRange([false, false]);

        // Add equality constraint
        g.Add(Xk_end - Xk);
        lbg.AddRange([0.0, 0.0]);
        ubg.AddRange([0.0, 0.0]);
      }
    }

    // Concatenate decision variables and constraint terms
    MX w_all = ca.vertcat(w);
    MX g_all = ca.vertcat(g);

    // Create an NLP solver
    var nlp_prob = new Dictionary<string, MX> { ["f"] = J, ["x"] = w_all, ["g"] = g_all };
    Function nlp_solver = ca.nlpsol("nlp_solver", "bonmin", nlp_prob, new Dictionary<string, object> { ["discrete"] = discrete });
    //Function nlp_solver = nlpsol("nlp_solver", "knitro", nlp_prob, new Dictionary<string, object> { ["discrete"] = discrete });
    //Function nlp_solver = nlpsol("nlp_solver", "ipopt", nlp_prob); // Solve relaxed problem

    // Plot the solution
    double[] tgrid = Enumerable.Range(0, N + 1).Select(k => T / N * k).ToArray();

    // Solve the NLP
    var sol = nlp_solver.call(new Dictionary<string, DM> {
      ["x0"] = ca.vertcat(w0), ["lbx"] = lbw.ToArray(), ["ubx"] = ubw.ToArray(),
      ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() });

    Console.WriteLine($"{{{string.Join(", ", nlp_solver.stats().Select(e => $"{e.Key}={e.Value}"))}}}");

    DM w1_opt = sol["x"];
    DM lam_w_opt = sol["lam_x"];
    DM lam_g_opt = sol["lam_g"];
    // plot: w1_opt[0::3], w1_opt[1::3] and (step) w1_opt[2::3] versus tgrid
  }
}
