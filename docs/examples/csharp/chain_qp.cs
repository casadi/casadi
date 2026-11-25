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
// -*- coding: utf-8 -*-
/*
We want to model a chain attached to two supports and hanging in between. Let us discretise
it with N mass points connected by N-1 springs. Each mass i has position (yi,zi), i=1,...,N.
The equilibrium point of the system minimises the potential energy.

The potential energy of each spring is
Vi=D_i/2 * ((y_i-y_{i+1})^2 + (z_i-z_{i+1})^2)

The gravitational potential energy of each mass is
Vg_i = m_i*g0*z_i

The total potential energy is thus given by:

Vchain(y,z) = 1/2*sum{i=1,...,N-1} D_i ((y_i-y_{i+1})^2+(z_i-z_{i+1})^2) + g0 * sum{i=1,...,N} m_i * z_i

where y=[y_1,...,y_N] and z=[z_1,...,z_N]

We wish to solve
minimize{y,z} Vchain(y, z)

Subject to the piecewise linear ground constraints:
z_i >= zin
z_i - 0.1*y_i >= 0.5
*/
using System;
using System.Collections.Generic;
using CasADi;

public static class chain_qp {
  public static void Main() {
    if (!ca.has_conic("qpoases")) { Console.WriteLine("SKIP: conic qpoases not available"); return; }

    // Constants
    int N = 40;
    double m_i = 40.0 / N;
    double D_i = 70.0 * N;
    double g0 = 9.81;
    //double zmin = -inf; // unbounded
    double zmin = 0.5; // ground

    // Objective function
    SX Vchain = 0;

    // Variables
    var x = new List<SX>();

    // Variable bounds
    var lbx = new List<double>();
    var ubx = new List<double>();

    // Constraints
    var g = new List<SX>();

    // Constraint bounds
    var lbg = new List<double>();
    var ubg = new List<double>();

    // Loop over all chain elements
    SX y_i = null, z_i = null, y_prev = null, z_prev = null;
    for (int i = 1; i <= N; i++) {
      // Previous point
      if (i > 1) {
        y_prev = y_i;
        z_prev = z_i;
      }

      // Create variables for the (y_i, z_i) coordinates
      y_i = SX.sym($"y_{i}");
      z_i = SX.sym($"z_{i}");

      // Add to the list of variables
      x.AddRange([y_i, z_i]);
      if (i == 1) {
        lbx.AddRange([-2.0, 1.0]);
        ubx.AddRange([-2.0, 1.0]);
      } else if (i == N) {
        lbx.AddRange([2.0, 1.0]);
        ubx.AddRange([2.0, 1.0]);
      } else {
        lbx.AddRange([-ca.inf, zmin]);
        ubx.AddRange([ca.inf, ca.inf]);
      }

      // Spring potential
      if (i > 1) {
        Vchain += D_i / 2 * (ca.pow(y_prev - y_i, 2) + ca.pow(z_prev - z_i, 2));
      }

      // Graviational potential
      Vchain += g0 * m_i * z_i;

      // Slanted ground constraints
      g.Add(z_i - 0.1 * y_i);
      lbg.Add(0.5);
      ubg.Add(ca.inf);
    }

    // Formulate QP
    var qp = new Dictionary<string, SX> { ["x"] = ca.vertcat(x), ["f"] = Vchain, ["g"] = ca.vertcat(g) };

    // Solve with IPOPT
    Function solver = ca.qpsol("solver", "qpoases", qp, new Dictionary<string, object> { ["sparse"] = true });
    //Function solver = qpsol("solver", "gurobi", qp);
    //Function solver = nlpsol("solver", "ipopt", qp);

    // Get the optimal solution
    var sol = solver.call(new Dictionary<string, DM> {
      ["lbx"] = lbx.ToArray(), ["ubx"] = ubx.ToArray(),
      ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() });
    DM x_opt = sol["x"];
    DM f_opt = sol["f"];
    Console.WriteLine($"f_opt =  {f_opt}");

    // Retrieve the result
    DM Y0 = x_opt[ca.Slice(0, x_opt.numel(), 2)];
    DM Z0 = x_opt[ca.Slice(1, x_opt.numel(), 2)];

    // Plot the result
    DM ys = ca.linspace(-2.0, 2.0, 100);
    DM zs = 0.5 + 0.1 * ys;
    // plot: chain (Y0, Z0) and ground line (ys, zs), 'hanging chain QP'
  }
}
