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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;

/*
Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
@author: Joel Andersson, KU Leuven 2013
*/
public class implicit_runge_kutta {
  public static void main(String[] args) {
    if (!ca.has_rootfinder("newton")) { System.out.println("SKIP: rootfinder newton not available"); return; }
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }

    // End time
    double tf = 10.0;

    // Dimensions
    int nx = 3;
    int np = 1;

    // Declare variables
    SX x = SX.sym("x", nx);  // state
    SX p = SX.sym("u", np);  // control

    // ODE right hand side function
    SX ode = ca.vertcat(List.of(ca.minus(1, x.get(1).times(x.get(1))).times(x.get(0)).minus(x.get(1)).plus(p),
        x.get(0),
        x.get(0).times(x.get(0)).plus(x.get(1).times(x.get(1))).plus(p.times(p))));
    Map<String, SX> dae = Map.of("x", x, "p", p, "ode", ode);
    Function f = new Function("f", List.of(x, p), List.of(ode));

    // Number of finite elements
    int n = 100;

    // Size of the finite elements
    double h = tf / n;

    // Degree of interpolating polynomial
    int d = 4;

    // Choose collocation points
    double[] tau_root = new double[d + 1];
    System.arraycopy(ca.collocation_points(d, "legendre"), 0, tau_root, 1, d);

    // Coefficients of the collocation equation
    double[][] C = new double[d + 1][d + 1];

    // Coefficients of the continuity equation
    double[] D = new double[d + 1];

    // Dimensionless time inside one control interval
    SX tau = SX.sym("tau");

    // For all collocation points
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      SX L = new SX(1);
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          L = L.times(tau.minus(tau_root[r]).div(tau_root[j] - tau_root[r]));
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      Function lfcn = new Function("lfcn", List.of(tau), List.of(L));
      D[j] = ((DM) lfcn.invoke(1.0)).to_double();

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      Function tfcn = new Function("tfcn", List.of(tau), List.of(ca.tangent(L, tau)));
      for (int r = 0; r < d + 1; r++) C[j][r] = ((DM) tfcn.invoke(tau_root[r])).to_double();
    }

    // Total number of variables for one finite element
    MX X0 = MX.sym("X0", nx);
    MX P = MX.sym("P", np);
    MX V = MX.sym("V", d * nx);

    // Get the state at each collocation point
    long[] offsets = new long[d + 1];
    for (int r = 0; r < d + 1; r++) offsets[r] = r * nx;
    List<MX> X = new ArrayList<>(List.of(X0));
    X.addAll(ca.vertsplit(V, offsets));

    // Get the collocation quations (that define V)
    List<MX> V_eq = new ArrayList<>();
    for (int j = 1; j < d + 1; j++) {
      // Expression for the state derivative at the collocation point
      MX xp_j = new MX(0);
      for (int r = 0; r < d + 1; r++) {
        xp_j = xp_j.plus(ca.times(C[r][j], X.get(r)));
      }

      // Append collocation equations
      MX f_j = f.call(List.of(X.get(j), P)).get(0);
      V_eq.add(ca.times(h, f_j).minus(xp_j));
    }

    // Concatenate constraints
    MX V_eq_all = ca.vertcat(V_eq);

    // Root-finding function, implicitly defines V as a function of X0 and P
    Function vfcn = new Function("vfcn", List.of(V, X0, P), List.of(V_eq_all));

    // Convert to SX to decrease overhead
    Function vfcn_sx = vfcn.expand();

    // Create a implicit function instance to solve the system of equations
    Function ifcn = ca.rootfinder("ifcn", "newton", vfcn_sx);
    V = ifcn.call(List.of(new MX(), X0, P)).get(0);
    X = new ArrayList<>();
    for (int r = 0; r < d + 1; r++) X.add(r == 0 ? X0 : V.get(new Slice((r - 1) * nx, r * nx)));

    // Get an expression for the state at the end of the finie element
    MX XF = new MX(0);
    for (int r = 0; r < d + 1; r++) {
      XF = XF.plus(ca.times(D[r], X.get(r)));
    }

    // Get the discrete time dynamics
    Function F = new Function("F", List.of(X0, P), List.of(XF));

    // Do this iteratively for all finite elements
    MX Xi = X0;
    for (int i = 0; i < n; i++) {
      Xi = F.call(List.of(Xi, P)).get(0);
    }

    // Fixed-step integrator
    Function irk_integrator = new Function("irk_integrator", Map.of("x0", X0, "p", P, "xf", Xi),
        ca.integrator_in(), ca.integrator_out());

    // Create a convensional integrator for reference
    Function ref_integrator = ca.integrator("ref_integrator", "cvodes", dae, 0, tf);

    // Test values
    double[] x0_val = {0, 1, 0};
    double p_val = 0.2;

    // Make sure that both integrators give consistent results
    for (Function Fi : List.of(irk_integrator, ref_integrator)) {
      System.out.println("-------");
      System.out.println("Testing " + Fi.name());
      System.out.println("-------");

      // Generate a new function that calculates forward and reverse directional derivatives
      Function dF = Fi.factory("dF", List.of("x0", "p", "fwd:x0", "fwd:p", "adj:xf"),
                                     List.of("xf", "fwd:xf", "adj:x0", "adj:p"));
      Map<String, Object> arg = new LinkedHashMap<>();

      // Pass arguments
      arg.put("x0", x0_val);
      arg.put("p", p_val);

      // Forward sensitivity analysis, first direction: seed p and x0[0]
      arg.put("fwd_x0", new double[] {1, 0, 0});
      arg.put("fwd_p", 1);

      // Adjoint sensitivity analysis, seed xf[2]
      arg.put("adj_xf", new double[] {0, 0, 1});

      // Integrate
      Map<String, GenericMatrix> res = dF.invoke(arg);

      // Get the nondifferentiated results
      System.out.println(String.format("%30s = %s", "xf", res.get("xf")));

      // Get the forward sensitivities
      System.out.println(String.format("%30s = %s", "d(xf)/d(p)+d(xf)/d(x0[0])", res.get("fwd_xf")));

      // Get the adjoint sensitivities
      System.out.println(String.format("%30s = %s", "d(xf[2])/d(x0)", res.get("adj_x0")));
      System.out.println(String.format("%30s = %s", "d(xf[2])/d(p)", res.get("adj_p")));
    }
  }
}
