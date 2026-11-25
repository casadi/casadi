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
// -*- coding: utf-8 -*-

using System;
using System.Collections.Generic;
using CasADi;

public static class sensitivity_analysis {
  // Python print('%50s' % label, a, x, b, y)
  static void print(string label, string a, object x, string b, object y) {
    Console.WriteLine($"{label,50} {a} {x} {b} {y}");
  }

  public static void Main() {
    foreach (string p in new[] { "cvodes", "rk", "idas", "collocation" }) {
      if (!ca.has_integrator(p)) { Console.WriteLine($"SKIP: integrator {p} not available"); return; }
    }
    if (!ca.has_rootfinder("kinsol")) { Console.WriteLine("SKIP: rootfinder kinsol not available"); return; }
    Console.WriteLine("Testing sensitivity analysis in CasADi");

    // All ODE and DAE integrators to be tested
    var DAE_integrators = new List<string> { "idas", "collocation" };
    var ODE_integrators = new List<string> { "cvodes", "rk", "idas", "collocation" };

    foreach (var Integrators in new[] { ODE_integrators, DAE_integrators }) {
      Dictionary<string, SX> dae;
      double tf;
      DM x0;
      double u0;
      if (Integrators == ODE_integrators) { // rocket example
        Console.WriteLine("******");
        Console.WriteLine("Testing ODE example");

        // Time
        SX t = SX.sym("t");

        // Parameter
        SX u = SX.sym("u");

        // Differential states
        SX s = SX.sym("s"), v = SX.sym("v"), m = SX.sym("m");
        SX x = ca.vertcat(s, v, m);

        // Constants
        double alpha = 0.05; // friction
        double beta = 0.1;   // fuel consumption rate

        // Differential equation
        SX ode = ca.vertcat(
          v,
          (u - alpha * v * v) / m,
          -beta * u * u);

        // Quadrature
        SX quad = ca.pow(v, 3) + ca.pow(3 - ca.sin(t) - u, 2);

        // DAE callback function
        dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = u, ["ode"] = ode, ["quad"] = quad };

        // Time length
        tf = 0.5;

        // Initial position
        x0 = new double[] { 0.0, 0.0, 1.0 };

        // Parameter
        u0 = 0.4;

      } else { // Simple DAE example
        Console.WriteLine("******");
        Console.WriteLine("Testing DAE example");

        // Differential state
        SX x = SX.sym("x");

        // Algebraic variable
        SX z = SX.sym("z");

        // Parameter
        SX u = SX.sym("u");

        // Differential equation
        SX ode = -x + 0.5 * x * x + u + 0.5 * z;

        // Algebraic constraint
        SX alg = z + ca.exp(z) - 1.0 + x;

        // Quadrature
        SX quad = x * x + 3.0 * u * u;

        // DAE callback function
        dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = u, ["ode"] = ode, ["alg"] = alg, ["quad"] = quad };

        // End time
        tf = 5.0;

        // Initial position
        x0 = 1.0;

        // Parameter
        u0 = 0.4;
      }

      // Integrator
      foreach (string MyIntegrator in Integrators) {
        Console.WriteLine("========");
        Console.WriteLine($"Integrator:  {MyIntegrator}");
        Console.WriteLine("========");

        // Integrator options
        var opts = new Dictionary<string, object>();
        if (MyIntegrator == "collocation") {
          opts["rootfinder"] = "kinsol";
        }

        // Integrator
        Function I = ca.integrator("I", MyIntegrator, dae, 0, tf, opts);

        // Integrate to get results
        var res = I.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0 });
        DM xf = res["xf"];
        DM qf = res["qf"];
        print("Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf);

        // Perturb solution to get a finite difference approximation
        double h = 0.001;
        res = I.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0 + h });
        DM fd_xf = (res["xf"] - xf) / h;
        DM fd_qf = (res["qf"] - qf) / h;
        print("Finite difference approximation:", "d(xf)/d(p) = ", fd_xf, ", d(qf)/d(p) = ", fd_qf);

        // Calculate once, forward
        Function I_fwd = I.factory("I_fwd", ["x0", "z0", "p", "fwd:p"], ["fwd:xf", "fwd:qf"]);
        res = I_fwd.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0, ["fwd_p"] = 1 });
        DM fwd_xf = res["fwd_xf"];
        DM fwd_qf = res["fwd_qf"];
        print("Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf);

        // Calculate once, adjoint
        Function I_adj = I.factory("I_adj", ["x0", "z0", "p", "adj:qf"], ["adj:x0", "adj:p"]);
        res = I_adj.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0, ["adj_qf"] = 1 });
        DM adj_x0 = res["adj_x0"];
        DM adj_p = res["adj_p"];
        print("Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p);

        // Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
        res = I_adj.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0 + h, ["adj_qf"] = 1 });
        DM fd_adj_x0 = (res["adj_x0"] - adj_x0) / h;
        DM fd_adj_p = (res["adj_p"] - adj_p) / h;
        print("FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fd_adj_x0, ", d2(qf)/d(p)d(p) = ", fd_adj_p);

        // Forward over adjoint to get the second order sensitivities
        Function I_foa = I_adj.factory("I_foa", ["x0", "z0", "p", "adj_qf", "fwd:p"], ["fwd:adj_x0", "fwd:adj_p"]);
        res = I_foa.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0, ["adj_qf"] = 1, ["fwd_p"] = 1 });
        DM fwd_adj_x0 = res["fwd_adj_x0"];
        DM fwd_adj_p = res["fwd_adj_p"];
        print("Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p);

        // Adjoint over adjoint to get the second order sensitivities
        Function I_aoa = I_adj.factory("I_aoa", ["x0", "z0", "p", "adj_qf", "adj:adj_p"], ["adj:x0", "adj:p"]);
        res = I_aoa.call(new Dictionary<string, DM> { ["x0"] = x0, ["p"] = u0, ["adj_qf"] = 1, ["adj_adj_p"] = 1 });
        adj_x0 = res["adj_x0"];
        adj_p = res["adj_p"];
        print("Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_x0, ", d2(qf)/d(p)d(p) = ", adj_p);
      }
    }
  }
}
