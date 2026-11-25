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

// An implementation of direct collocation
// Joris Gillis, 2018
public static class direct_collocation_opti {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Degree of interpolating polynomial
    int d = 3;

    // Get collocation points
    List<double> tau = ca.collocation_points(d, "legendre");

    // Collocation linear maps
    var (C, D, B) = ca.collocation_coeff(tau);

    // Time horizon
    double T = 10;

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
    var f = ca.Function("f", [x, u], [xdot, L]);

    // Control discretization
    int N = 20; // number of control intervals
    double h = T / N;

    // Start with an empty NLP

    var opti = ca.Opti();
    MX J = 0;

    // "Lift" initial conditions
    MX Xk = opti.variable(2);
    opti.subject_to(ca.eq(Xk, new double[] { 0, 1 }));
    opti.set_initial(Xk, new double[] { 0, 1 });

    // Collect all states/controls
    var Xs = new List<MX> { Xk };
    var Us = new List<MX>();

    // Formulate the NLP
    for (int k = 0; k <= N - 1; k++) {
      // New NLP variable for the control
      MX Uk = opti.variable();
      Us.Add(Uk);
      opti.subject_to(-1 <= (Uk <= 1));
      opti.set_initial(Uk, 0);

      // Decision variables for helper states at each collocation point
      MX Xc = opti.variable(2, d);
      opti.subject_to(-0.25 <= Xc[0, ..]);
      opti.set_initial(Xc, ca.repmat(new double[] { 0, 0 }, 1, d));

      // Evaluate ODE right-hand-side at all helper states
      List<MX> odequad = f.call([Xc, Uk]);
      MX ode = odequad[0], quad = odequad[1];

      // Add contribution to quadrature function
      J += ca.mtimes(quad, B) * h;

      // Get interpolating points of collocation polynomial
      MX Z = ca.horzcat(Xk, Xc);

      // Get slope of interpolating polynomial (normalized)
      MX Pidot = ca.mtimes(Z, C);
      // Match with ODE right-hand-side
      opti.subject_to(ca.eq(Pidot, h * ode));

      // State at end of collocation interval
      MX Xk_end = ca.mtimes(Z, D);

      // New decision variable for state at end of interval
      Xk = opti.variable(2);
      Xs.Add(Xk);
      opti.subject_to(-0.25 <= Xk[0]);
      opti.set_initial(Xk, new double[] { 0, 0 });

      // Continuity constraints
      opti.subject_to(ca.eq(Xk_end, Xk));
    }

    MX Xs_ = ca.horzcat(Xs);
    MX Us_ = ca.horzcat(Us);

    opti.minimize(J);

    opti.solver("ipopt");

    OptiSol sol = opti.solve();

    DM x_opt = sol.value(Xs_);
    DM u_opt = sol.value(Us_);

    // Plot the solution
    DM tgrid = ca.linspace(0.0, T, N + 1);
    // plot: x_opt(1,:) dashed and x_opt(2,:) solid against tgrid, stairs of [u_opt nan] dash-dotted;
    //       xlabel t, legend x1, x2, u
  }
}
