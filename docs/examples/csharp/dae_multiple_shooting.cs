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
//
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using CasADi;

/*
Solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE)

minimize     integral_{t=0}^{10} x0^2 + x1^2 + u^2  dt
x0,x1,z,u

subject to   dot(x0) == z*x0-x1+u     \
             dot(x1) == x0             }  for 0 <= t <= 10
                   0 == x1^2 + z - 1  /
             x0(t=0) == 0
             x1(t=0) == 1
             x0(t=10) == 0
             x1(t=10) == 0
             -0.75 <= u <= 1  for 0 <= t <= 10

The method used is direct multiple shooting.

Joel Andersson, 2015
*/
public static class dae_multiple_shooting {
  public static void Main() {
    if (!ca.has_integrator("idas")) { Console.WriteLine("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    if (!hasLibrary("hsl")) { Console.WriteLine("SKIP: libhsl (ipopt linear_solver ma27) not available"); return; }

    // Declare variables
    SX x0 = SX.sym("x0");
    SX x1 = SX.sym("x1");
    SX x = ca.vertcat(x0, x1); // Differential states
    SX z = SX.sym("z");       // Algebraic variable
    SX u = SX.sym("u");       // Control

    // Differential equation
    SX f_x = ca.vertcat(z * x0 - x1 + u, x0);

    // Algebraic equation
    SX f_z = ca.pow(x1, 2) + z - 1;

    // Lagrange cost term (quadrature)
    SX f_q = ca.pow(x0, 2) + ca.pow(x1, 2) + ca.pow(u, 2);

    // Create an integrator
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = u, ["ode"] = f_x, ["alg"] = f_z, ["quad"] = f_q };
    // interval length 0.5s
    Function I = ca.integrator("I", "idas", dae, 0, 0.5);

    // Number of intervals
    int nk = 20;

    // Start with an empty NLP
    var w = new List<MX>();         // List of variables
    var lbw = new List<double>();   // Lower bounds on w
    var ubw = new List<double>();   // Upper bounds on w
    var G = new List<MX>();         // Constraints
    MX J = 0;                       // Cost function

    // Initial conditions
    MX Xk = MX.sym("X0", 2);
    w.Add(Xk);
    lbw.AddRange([0.0, 1.0]);
    ubw.AddRange([0.0, 1.0]);

    // Loop over all intervals
    for (int k = 0; k < nk; k++) {
      // Local control
      MX Uk = MX.sym($"U{k}");
      w.Add(Uk);
      lbw.Add(-0.75);
      ubw.Add(1.00);

      // Call integrator function
      var Ik = I.call(new Dictionary<string, MX> { ["x0"] = Xk, ["p"] = Uk });
      Xk = Ik["xf"];
      J += Ik["qf"]; // Sum quadratures

      // "Lift" the variable
      MX X_prev = Xk;
      Xk = MX.sym($"X{k + 1}", 2);
      w.Add(Xk);
      lbw.AddRange([-ca.inf, -ca.inf]);
      ubw.AddRange([ca.inf, ca.inf]);
      G.Add(X_prev - Xk);
    }

    // Allocate an NLP solver
    var nlp = new Dictionary<string, MX> { ["x"] = ca.vertcat(w), ["f"] = J, ["g"] = ca.vertcat(G) };
    var opts = new Dictionary<string, object> { ["ipopt.linear_solver"] = "ma27" };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Pass bounds, initial guess and solve NLP
    var sol = solver.call(new Dictionary<string, DM> {
      ["lbx"] = lbw.ToArray(), // Lower variable bound
      ["ubx"] = ubw.ToArray(), // Upper variable bound
      ["lbg"] = 0.0,           // Lower constraint bound
      ["ubg"] = 0.0,           // Upper constraint bound
      ["x0"] = 0.0 });         // Initial guess

    // Plot the results
    DM x_opt = sol["x"];
    DM tgrid = ca.linspace(0.0, 10.0, nk + 1);
    DM x0_opt = x_opt[ca.Slice(0, x_opt.numel(), 3)];
    DM x1_opt = x_opt[ca.Slice(1, x_opt.numel(), 3)];
    DM u_opt = x_opt[ca.Slice(2, x_opt.numel(), 3)];
    // plot: x0, x1 (vs tgrid) and u trajectories, 'Van der Pol optimization - multiple shooting'
  }

  static bool hasLibrary(string name) {
    if (!NativeLibrary.TryLoad(name, typeof(dae_multiple_shooting).Assembly, null, out IntPtr handle)) return false;
    NativeLibrary.Free(handle);
    return true;
  }
}
