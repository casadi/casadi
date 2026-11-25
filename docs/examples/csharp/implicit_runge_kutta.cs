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
using System.Linq;
using CasADi;

/*
Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
@author: Joel Andersson, KU Leuven 2013
*/
public static class implicit_runge_kutta {
  public static void Main() {
    if (!ca.has_rootfinder("newton")) { Console.WriteLine("SKIP: rootfinder newton not available"); return; }
    if (!ca.has_integrator("cvodes")) { Console.WriteLine("SKIP: integrator cvodes not available"); return; }

    // End time
    double tf = 10.0;

    // Dimensions
    int nx = 3;
    int np = 1;

    // Declare variables
    SX x = SX.sym("x", nx);  // state
    SX p = SX.sym("u", np);  // control

    // ODE right hand side function
    SX ode = ca.vertcat((1 - x[1] * x[1]) * x[0] - x[1] + p,
                     x[0],
                     x[0] * x[0] + x[1] * x[1] + p * p);
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ode };
    var f = ca.Function("f", [x, p], [ode]);

    // Number of finite elements
    int n = 100;

    // Size of the finite elements
    double h = tf / n;

    // Degree of interpolating polynomial
    int d = 4;

    // Choose collocation points
    double[] tau_root = ca.collocation_points(d, "legendre").Prepend(0.0).ToArray();

    // Coefficients of the collocation equation
    var C = new double[d + 1, d + 1];

    // Coefficients of the continuity equation
    var D = new double[d + 1];

    // Dimensionless time inside one control interval
    SX tau = SX.sym("tau");

    // For all collocation points
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      SX L = 1;
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          L *= (tau - tau_root[r]) / (tau_root[j] - tau_root[r]);
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      var lfcn = ca.Function("lfcn", [tau], [L]);
      D[j] = (double)lfcn.Call(1.0);

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      var tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)]);
      for (int r = 0; r < d + 1; r++) C[j, r] = (double)tfcn.Call(tau_root[r]);
    }

    // Total number of variables for one finite element
    MX X0 = MX.sym("X0", nx);
    MX P = MX.sym("P", np);
    MX V = MX.sym("V", d * nx);

    // Get the state at each collocation point
    var X = new List<MX> { X0 };
    X.AddRange(ca.vertsplit(V, Enumerable.Range(0, d + 1).Select(r => (long)(r * nx))));

    // Get the collocation quations (that define V)
    var V_eq = new List<MX>();
    for (int j = 1; j < d + 1; j++) {
      // Expression for the state derivative at the collocation point
      MX xp_j = 0;
      for (int r = 0; r < d + 1; r++) {
        xp_j += C[r, j] * X[r];
      }

      // Append collocation equations
      MX f_j = f.Call(X[j], P);
      V_eq.Add(h * f_j - xp_j);
    }

    // Concatenate constraints
    MX V_eq_all = ca.vertcat(V_eq);

    // Root-finding function, implicitly defines V as a function of X0 and P
    var vfcn = ca.Function("vfcn", [V, X0, P], [V_eq_all]);

    // Convert to SX to decrease overhead
    Function vfcn_sx = vfcn.expand();

    // Create a implicit function instance to solve the system of equations
    Function ifcn = ca.rootfinder("ifcn", "newton", vfcn_sx);
    V = ifcn.Call(ca.MX(), X0, P);
    X = new List<MX>();
    for (int r = 0; r < d + 1; r++) X.Add(r == 0 ? X0 : V[((r - 1) * nx)..(r * nx)]);

    // Get an expression for the state at the end of the finie element
    MX XF = 0;
    for (int r = 0; r < d + 1; r++) {
      XF += D[r] * X[r];
    }

    // Get the discrete time dynamics
    var F = ca.Function("F", [X0, P], [XF]);

    // Do this iteratively for all finite elements
    MX Xi = X0;
    for (int i = 0; i < n; i++) {
      Xi = F.Call(Xi, P);
    }

    // Fixed-step integrator
    var irk_integrator = ca.Function("irk_integrator", new Dictionary<string, MX> { ["x0"] = X0, ["p"] = P, ["xf"] = Xi },
        ca.integrator_in(), ca.integrator_out());

    // Create a convensional integrator for reference
    Function ref_integrator = ca.integrator("ref_integrator", "cvodes", dae, 0, tf);

    // Test values
    double[] x0_val = { 0, 1, 0 };
    double p_val = 0.2;

    // Make sure that both integrators give consistent results
    foreach (Function Fi in new[] { irk_integrator, ref_integrator }) {
      Console.WriteLine("-------");
      Console.WriteLine($"Testing {Fi.name()}");
      Console.WriteLine("-------");

      // Generate a new function that calculates forward and reverse directional derivatives
      Function dF = Fi.factory("dF", ["x0", "p", "fwd:x0", "fwd:p", "adj:xf"],
                                     ["xf", "fwd:xf", "adj:x0", "adj:p"]);
      var arg = new Dictionary<string, DM>();

      // Pass arguments
      arg["x0"] = x0_val;
      arg["p"] = p_val;

      // Forward sensitivity analysis, first direction: seed p and x0[0]
      arg["fwd_x0"] = new double[] { 1, 0, 0 };
      arg["fwd_p"] = 1;

      // Adjoint sensitivity analysis, seed xf[2]
      arg["adj_xf"] = new double[] { 0, 0, 1 };

      // Integrate
      var res = dF.call(arg);

      // Get the nondifferentiated results
      Console.WriteLine($"{"xf",30} = {res["xf"]}");

      // Get the forward sensitivities
      Console.WriteLine($"{"d(xf)/d(p)+d(xf)/d(x0[0])",30} = {res["fwd_xf"]}");

      // Get the adjoint sensitivities
      Console.WriteLine($"{"d(xf[2])/d(x0)",30} = {res["adj_x0"]}");
      Console.WriteLine($"{"d(xf[2])/d(p)",30} = {res["adj_p"]}");
    }
  }
}
