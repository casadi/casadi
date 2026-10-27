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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.ca;

public class biegler_10_1 {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Excercise 1, chapter 10 from Larry Biegler's book
    System.out.println("program started");

    // Test with different number of elements
    for (int N = 1; N < 11; N++) {
      System.out.println("N =  " + N);

      // Degree of interpolating polynomial
      int K = 2;

      // Legrandre roots
      double[] tau_root = {0., 0.211325, 0.788675};

      // Radau roots (K=3)
      //double[] tau_root = {0, 0.155051, 0.644949, 1};

      // Time
      SX t = SX.sym("t");

      // Differential equation
      SX z = SX.sym("z");
      Function F = new Function("dz_dt", List.of(z), List.of(z.times(z).minus(z.times(2)).plus(1)));

      double z0 = -3;

      // Analytic solution
      Function z_analytic = new Function("z_analytic", List.of(t), List.of(t.times(4).minus(3).div(t.times(3).plus(1))));

      // Collocation point
      SX tau = SX.sym("tau");

      // Step size
      double h = 1.0 / N;

      // Get the coefficients of the continuity and collocation equations
      DM D = DM.zeros(K + 1);
      DM C = DM.zeros(K + 1, K + 1);
      for (int j = 0; j < K + 1; j++) {
        // Lagrange polynomial
        SX L = new SX(1);
        for (int k = 0; k < K + 1; k++) {
          if (k != j) {
            L = L.times(tau.minus(tau_root[k]).div(tau_root[j] - tau_root[k]));
          }
        }

        // Evaluate at end for coefficients of continuity equation
        Function lfcn = new Function("lfcn", List.of(tau), List.of(L));
        D.set(j, (DM) lfcn.invoke(1.));

        // Differentiate and evaluate at collocation points
        Function tfcn = new Function("tfcn", List.of(tau), List.of(ca.tangent(L, tau)));
        for (int k = 0; k < K + 1; k++) C.set(j, k, (DM) tfcn.invoke(tau_root[k]));
      }
      System.out.println("C =  " + C);
      System.out.println("D =  " + D);

      // Collocated states
      SX Z = SX.sym("Z", N, K + 1);

      // Construct the NLP
      SX x = ca.vec(Z.T());
      List<SX> g = new ArrayList<>();
      for (int i = 0; i < N; i++) {
        for (int k = 1; k < K + 1; k++) {
          // Add collocation equations to NLP
          SX rhs = new SX(0);
          for (int j = 0; j < K + 1; j++) {
            rhs = rhs.plus(Z.get(i, j).times(C.get(j, k)));
          }
          SX FF = F.call(List.of(Z.get(i, k))).get(0);
          g.add(ca.times(h, FF).minus(rhs));
        }

        // Add continuity equation to NLP
        SX rhs = new SX(0);
        for (int j = 0; j < K + 1; j++) {
          rhs = rhs.plus(ca.times(D.get(j), Z.get(i, j)));
        }

        if (i < N - 1) {
          g.add(Z.get(i + 1, 0).minus(rhs));
        }
      }

      SX gg = ca.vertcat(g);

      System.out.println("g =  " + gg);

      // NLP
      Map<String, SX> nlp = Map.of("x", x, "f", x.get(0).power(2), "g", gg);

      //# ----
      //# SOLVE THE NLP
      //# ----

      // NLP solver options
      Map<String, Object> opts = Map.of("ipopt.tol", 1e-10);

      // Allocate an NLP solver and buffer
      Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
      Map<String, Object> arg = new LinkedHashMap<>();

      // Initial condition
      arg.put("x0", new double[(int) x.nnz()]);

      // Bounds on x
      double[] lbx = new double[(int) x.nnz()];
      double[] ubx = new double[(int) x.nnz()];
      Arrays.fill(lbx, -100);
      Arrays.fill(ubx, 100);
      lbx[0] = ubx[0] = z0;
      arg.put("lbx", lbx);
      arg.put("ubx", ubx);

      // Bounds on the constraints
      arg.put("lbg", 0);
      arg.put("ubg", 0);

      // Solve the problem
      Map<String, GenericMatrix> res = solver.invoke(arg);

      //# Print the time points
      double[] t_opt = new double[N * (K + 1)];
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < K + 1; j++) {
          t_opt[j + (K + 1) * i] = h * (i + tau_root[j]);
        }
      }

      System.out.println("time points:  " + Arrays.toString(t_opt));

      // Print the optimal cost
      System.out.println("optimal cost:  " + ((DM) res.get("f")).to_double());

      // Print the optimal solution
      double[] xopt = ((DM) res.get("x")).nonzeros();
      System.out.println("optimal solution:  " + Arrays.toString(xopt));

      // plot: xopt versus t_opt
    }

    // show the plots
  }
}
