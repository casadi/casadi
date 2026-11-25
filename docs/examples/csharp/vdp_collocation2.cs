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

public static class vdp_collocation2 {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { Console.WriteLine("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }

    int nk = 20;      // Control discretization
    double tf = 10.0; // End time

    // Declare variables (use scalar graph)
    SX t = SX.sym("t"); // time
    SX u = SX.sym("u"); // control

    // states = [x (vdp oscillator states, 2); L (helper state: Langrange integrand)]
    SX states = SX.sym("states", 3);

    // Create a structure for the right hand side
    SX x = states[0..2];
    SX rhs = ca.vertcat(
        ca.vertcat((1 - x[1] * x[1]) * x[0] - x[1] + u, x[0]),
        x[0] * x[0] + x[1] * x[1] + u * u);

    // ODE right hand side function
    var f = ca.Function("f", [t, states, u], [rhs]);

    // Objective function (meyer term)
    var m = ca.Function("m", [t, states, u], [states[2]]);

    // Control bounds
    double u_min = -0.75;
    double u_max = 1.0;
    double u_init = 0.0;

    // State bounds and initial guess
    double[] x_min = { -ca.inf, -ca.inf, -ca.inf };
    double[] x_max = {  ca.inf,  ca.inf,  ca.inf };
    double[] xi_min = { 0.0, 1.0, 0.0 };
    double[] xi_max = { 0.0, 1.0, 0.0 };
    double[] xf_min = { 0.0, 0.0, -ca.inf };
    double[] xf_max = { 0.0, 0.0,  ca.inf };
    double[] x_init = { 0.0, 0.0, 0.0 };

    // Dimensions
    int nx = 3;
    int nu = 1;

    // Choose collocation points
    double[] tau_root = ca.collocation_points(3, "radau").Prepend(0).ToArray();

    // Degree of interpolating polynomial
    int d = tau_root.Length - 1;

    // Size of the finite elements
    double h = tf / nk;

    // Coefficients of the collocation equation
    var C = new double[d + 1, d + 1];

    // Coefficients of the continuity equation
    var D = new double[d + 1];

    // Dimensionless time inside one control interval
    SX tau = SX.sym("tau");

    // All collocation time points
    var T = new double[nk, d + 1];
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < d + 1; j++)
        T[k, j] = h * (k + tau_root[j]);

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
      for (int r = 0; r < d + 1; r++) {
        C[j, r] = (double)tfcn.Call(tau_root[r]);
      }
    }

    // Structure holding NLP variables: V = [X[0][0..d], U[0], X[1][0..d], U[1], ..., X[nk][0..d]]
    int blk = (d + 1) * nx + nu;
    int NV = nk * blk + (d + 1) * nx;
    MX V = MX.sym("V", NV);
    var iX = new int[nk + 1, d + 1];
    var iU = new int[nk];
    for (int k = 0; k < nk + 1; k++)
      for (int j = 0; j < d + 1; j++) iX[k, j] = k * blk + j * nx;
    for (int k = 0; k < nk; k++) iU[k] = k * blk + (d + 1) * nx;

    // V["X",k,j]
    MX VX(int k, int j) => V[iX[k, j]..(iX[k, j] + nx)];

    var vars_lb = new double[NV];
    var vars_ub = new double[NV];
    var vars_init = new double[NV];

    // Set states and its bounds
    for (int k = 0; k < nk + 1; k++)
      for (int j = 0; j < d + 1; j++)
        for (int i = 0; i < nx; i++) {
          vars_init[iX[k, j] + i] = x_init[i];
          vars_lb[iX[k, j] + i] = x_min[i];
          vars_ub[iX[k, j] + i] = x_max[i];
        }

    // Set controls and its bounds
    for (int k = 0; k < nk; k++) {
      vars_init[iU[k]] = u_init;
      vars_lb[iU[k]] = u_min;
      vars_ub[iU[k]] = u_max;
    }

    // State at initial time
    for (int i = 0; i < nx; i++) {
      vars_lb[iX[0, 0] + i] = xi_min[i];
      vars_ub[iX[0, 0] + i] = xi_max[i];
    }

    // State at end time
    for (int i = 0; i < nx; i++) {
      vars_lb[iX[nk, 0] + i] = xf_min[i];
      vars_ub[iX[nk, 0] + i] = xf_max[i];
    }

    // Constraint function for the NLP
    var g = new List<MX>();
    var lbg = new List<double>();
    var ubg = new List<double>();

    // For all finite elements
    for (int k = 0; k < nk; k++) {

      // For all collocation points
      for (int j = 1; j < d + 1; j++) {

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = 0;
        for (int r = 0; r < d + 1; r++) {
          xp_jk += C[r, j] * VX(k, r);
        }

        // Add collocation equations to the NLP
        MX fk = f.Call(T[k, j], VX(k, j), V[iU[k]]);
        g.Add(h * fk - xp_jk);
        for (int i = 0; i < nx; i++) { lbg.Add(0.0); ubg.Add(0.0); } // equality constraints
      }

      // Get an expression for the state at the end of the finite element
      MX xf_k = 0;
      for (int r = 0; r < d + 1; r++) {
        xf_k += D[r] * VX(k, r);
      }

      // Add continuity equation to NLP
      g.Add(VX(k + 1, 0) - xf_k);
      for (int i = 0; i < nx; i++) { lbg.Add(0.0); ubg.Add(0.0); }
    }

    // Concatenate constraints
    MX gc = ca.vertcat(g);

    // Objective function
    MX fobj = m.Call(T[nk - 1, d], VX(nk, 0), V[iU[nk - 1]]);

    // NLP
    var nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = fobj, ["g"] = gc };

    // ----
    // SOLVE THE NLP
    // ----

    // Set options
    var opts = new Dictionary<string, object>();
    opts["expand"] = true;
    // opts["ipopt.max_iter"] = 4L;
    opts["ipopt.linear_solver"] = "ma27";

    // Allocate an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    var arg = new Dictionary<string, DM>();

    // Initial condition
    arg["x0"] = vars_init;

    // Bounds on x
    arg["lbx"] = vars_lb;
    arg["ubx"] = vars_ub;

    // Bounds on g
    arg["lbg"] = lbg.ToArray();
    arg["ubg"] = ubg.ToArray();

    // Solve the problem
    var res = solver.call(arg);

    // Print the optimal cost
    Console.WriteLine($"optimal cost:  {(double)res["f"]}");

    // Retrieve the solution
    List<double> opt = res["x"].nonzeros();

    // Get values at the beginning of each finite element
    double[] x0_opt = new double[nk + 1], x1_opt = new double[nk + 1], x2_opt = new double[nk + 1];
    for (int k = 0; k < nk + 1; k++) {
      x0_opt[k] = opt[iX[k, 0]];
      x1_opt[k] = opt[iX[k, 0] + 1];
      x2_opt[k] = opt[iX[k, 0] + 2];
    }

    var u_opt = new double[nk];
    for (int k = 0; k < nk; k++) u_opt[k] = opt[iU[k]];

    DM tgrid = ca.linspace(0, tf, nk + 1);
    DM tgrid_u = ca.linspace(0, tf, nk);

    // plot: x0_opt, x1_opt over tgrid and u_opt (step) over tgrid_u, "Van der Pol optimization"
  }

  static bool hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", new Dictionary<string, SX> { ["x"] = z, ["f"] = z * z },
        new Dictionary<string, object> { ["ipopt.linear_solver"] = "ma27", ["ipopt.print_level"] = 0L, ["print_time"] = false });
    s.call(new Dictionary<string, DM> { ["x0"] = 1 });
    return (bool)s.stats()["success"];
  }
}
