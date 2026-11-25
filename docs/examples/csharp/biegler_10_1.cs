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
using System;
using System.Collections.Generic;
using System.Linq;
using CasADi;

public static class biegler_10_1 {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Excercise 1, chapter 10 from Larry Biegler's book
    Console.WriteLine("program started");

    // Test with different number of elements
    for (int N = 1; N < 11; N++) {
      Console.WriteLine($"N =  {N}");

      // Degree of interpolating polynomial
      int K = 2;

      // Legrandre roots
      double[] tau_root = { 0.0, 0.211325, 0.788675 };

      // Radau roots (K=3)
      //double[] tau_root = { 0, 0.155051, 0.644949, 1 };

      // Time
      SX t = SX.sym("t");

      // Differential equation
      SX z = SX.sym("z");
      var F = ca.Function("dz_dt", [z], [z * z - 2 * z + 1]);

      double z0 = -3;

      // Analytic solution
      var z_analytic = ca.Function("z_analytic", [t], [(4 * t - 3) / (3 * t + 1)]);

      // Collocation point
      SX tau = SX.sym("tau");

      // Step size
      double h = 1.0 / N;

      // Get the coefficients of the continuity and collocation equations
      DM D = DM.zeros(K + 1);
      DM C = DM.zeros(K + 1, K + 1);
      for (int j = 0; j < K + 1; j++) {
        // Lagrange polynomial
        SX L = 1;
        for (int k = 0; k < K + 1; k++) {
          if (k != j) {
            L *= (tau - tau_root[k]) / (tau_root[j] - tau_root[k]);
          }
        }

        // Evaluate at end for coefficients of continuity equation
        var lfcn = ca.Function("lfcn", [tau], [L]);
        D[j] = lfcn.Call(1.0);

        // Differentiate and evaluate at collocation points
        var tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)]);
        for (int k = 0; k < K + 1; k++) C[j, k] = tfcn.Call(tau_root[k]);
      }
      Console.WriteLine($"C =  {C}");
      Console.WriteLine($"D =  {D}");

      // Collocated states
      SX Z = SX.sym("Z", N, K + 1);

      // Construct the NLP
      SX x = ca.vec(Z.T());
      var g = new List<SX>();
      for (int i = 0; i < N; i++) {
        for (int k = 1; k < K + 1; k++) {
          // Add collocation equations to NLP
          SX rhs = 0;
          for (int j = 0; j < K + 1; j++) {
            rhs += Z[i, j] * C[j, k];
          }
          SX FF = F.Call(Z[i, k]);
          g.Add(h * FF - rhs);
        }

        // Add continuity equation to NLP
        {
          SX rhs = 0;
          for (int j = 0; j < K + 1; j++) {
            rhs += D[j] * Z[i, j];
          }

          if (i < N - 1) {
            g.Add(Z[i + 1, 0] - rhs);
          }
        }
      }

      SX gg = ca.vertcat(g);

      Console.WriteLine($"g =  {gg}");

      // NLP
      var nlp = new Dictionary<string, SX> { ["x"] = x, ["f"] = ca.pow(x[0], 2), ["g"] = gg };

      //# ----
      //# SOLVE THE NLP
      //# ----

      // NLP solver options
      var opts = new Dictionary<string, object> { ["ipopt.tol"] = 1e-10 };

      // Allocate an NLP solver and buffer
      Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
      var arg = new Dictionary<string, DM>();

      // Initial condition
      arg["x0"] = new double[x.nnz()];

      // Bounds on x
      double[] lbx = Enumerable.Repeat(-100.0, (int)x.nnz()).ToArray();
      double[] ubx = Enumerable.Repeat(100.0, (int)x.nnz()).ToArray();
      lbx[0] = ubx[0] = z0;
      arg["lbx"] = lbx;
      arg["ubx"] = ubx;

      // Bounds on the constraints
      arg["lbg"] = 0;
      arg["ubg"] = 0;

      // Solve the problem
      var res = solver.call(arg);

      //# Print the time points
      double[] t_opt = new double[N * (K + 1)];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < K + 1; j++) {
          t_opt[j + (K + 1) * i] = h * (i + tau_root[j]);
        }
      }

      Console.WriteLine($"time points:  [{string.Join(", ", t_opt)}]");

      // Print the optimal cost
      Console.WriteLine($"optimal cost:  {(double)res["f"]}");

      // Print the optimal solution
      List<double> xopt = res["x"].nonzeros();
      Console.WriteLine($"optimal solution:  [{string.Join(", ", xopt)}]");

      // plot: xopt versus t_opt
    }

    // show the plots
  }
}
