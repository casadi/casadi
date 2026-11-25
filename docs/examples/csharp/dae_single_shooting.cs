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

/*
This example mainly intended for CasADi presentations. It contains a compact
implementation of a direct single shooting method for DAEs using a minimal
number of CasADi concepts.

It solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE):

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

Note that other methods such as direct collocation or direct multiple shooting
are usually preferably to the direct single shooting method in practise.

Joel Andersson, 2012-2015
*/

public static class dae_single_shooting {
  public static void Main() {
    if (!ca.has_integrator("idas")) { Console.WriteLine("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x", 2); // Differential states
    SX z = SX.sym("z");    // Algebraic variable
    SX u = SX.sym("u");    // Control

    // Differential equation
    SX f_x = ca.vertcat(z * x[0] - x[1] + u, x[0]);

    // Algebraic equation
    SX f_z = ca.pow(x[1], 2) + z - 1;

    // Lagrange cost term (quadrature)
    SX f_q = ca.pow(x[0], 2) + ca.pow(x[1], 2) + ca.pow(u, 2);

    // Create an integrator
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = u, ["ode"] = f_x, ["alg"] = f_z, ["quad"] = f_q };
    // interval length 0.5s
    Function I = ca.integrator("I", "idas", dae, 0, 0.5);

    // All controls
    MX U = MX.sym("U", 20);

    // Construct graph of integrator calls
    MX X = new double[] { 0, 1 };
    MX J = 0;
    for (int k = 0; k < 20; k++) {
      var Ik = I.call(new Dictionary<string, MX> { ["x0"] = X, ["p"] = U[k] });
      X = Ik["xf"];
      J += Ik["qf"];   // Sum up quadratures
    }

    // Allocate an NLP solver
    var nlp = new Dictionary<string, MX> { ["x"] = U, ["f"] = J, ["g"] = X };
    var opts = new Dictionary<string, object> { ["ipopt.linear_solver"] = "ma27" };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Pass bounds, initial guess and solve NLP
    var sol = solver.call(new Dictionary<string, DM> {
      ["lbx"] = -0.75, // Lower variable bound
      ["ubx"] = 1.0,   // Upper variable bound
      ["lbg"] = 0.0,   // Lower constraint bound
      ["ubg"] = 0.0,   // Upper constraint bound
      ["x0"] = 0.0 }); // Initial guess
    Console.WriteLine($"{{{string.Join(", ", sol.Select(e => $"{e.Key}={e.Value}"))}}}");
  }
}
