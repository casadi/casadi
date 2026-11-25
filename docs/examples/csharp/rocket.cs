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
using CasADi;

public static class rocket {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Control
    MX u = MX.sym("u");

    // State
    MX x = MX.sym("x", 3);
    MX s = x[0]; // position
    MX v = x[1]; // speed
    MX m = x[2]; // mass

    // ODE right hand side
    MX sdot = v;
    MX vdot = (u - 0.05 * v * v) / m;
    MX mdot = -0.1 * u * u;
    MX xdot = ca.vertcat(sdot, vdot, mdot);

    // ODE right hand side function
    var f = ca.Function("f", [x, u], [xdot]);

    // Integrate with Explicit Euler over 0.2 seconds
    double dt = 0.01;  // Time step
    MX xj = x;
    for (int j = 0; j < 20; j++) {
      MX fj = f.Call(xj, u);
      xj = xj + dt * fj;
    }

    // Discrete time dynamics function
    var F = ca.Function("F", [x, u], [xj]);

    // Number of control segments
    int nu = 50;

    // Control for all segments
    MX U = MX.sym("U", nu);

    // Initial conditions
    MX X0 = new double[] { 0, 0, 1 };

    // Integrate over all intervals
    MX X = X0;
    for (int k = 0; k < nu; k++) {
      X = F.Call(X, U[k]);
    }

    // Objective function and constraints
    MX J = ca.mtimes(U.T(), U); // u'*u in Matlab
    MX G = X[0..2];          // x(1:2) in Matlab

    // NLP
    var nlp = new Dictionary<string, MX> { ["x"] = U, ["f"] = J, ["g"] = G };

    // Allocate an NLP solver
    var opts = new Dictionary<string, object> { ["ipopt.tol"] = 1e-10, ["expand"] = true };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    var arg = new Dictionary<string, DM>();

    // Bounds on u and initial condition
    arg["lbx"] = -0.5;
    arg["ubx"] = 0.5;
    arg["x0"] = 0.4;

    // Bounds on g
    arg["lbg"] = new double[] { 10, 0 };
    arg["ubg"] = new double[] { 10, 0 };

    // Solve the problem
    var res = solver.call(arg);

    // Get the solution
    // plot: res["x"] and res["lam_x"]
  }
}
