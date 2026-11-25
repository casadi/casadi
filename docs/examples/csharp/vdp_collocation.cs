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

public static class vdp_collocation {
  // Polynomials as coefficient arrays, highest power first (numpy.poly1d convention)
  static double[] polymul(double[] a, double[] b) {
    var c = new double[a.Length + b.Length - 1];
    for (int i = 0; i < a.Length; i++)
      for (int j = 0; j < b.Length; j++) c[i + j] += a[i] * b[j];
    return c;
  }

  static double polyval(double[] p, double x) {
    double r = 0;
    foreach (double c in p) r = r * x + c;
    return r;
  }

  static double[] polyder(double[] p) {
    int n = p.Length - 1;
    var r = new double[Math.Max(n, 1)];
    for (int i = 0; i < n; i++) r[i] = p[i] * (n - i);
    return r;
  }

  static double[] polyint(double[] p) {
    int n = p.Length;
    var r = new double[n + 1];
    for (int i = 0; i < n; i++) r[i] = p[i] / (n - i);
    return r;
  }

  static bool hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", new Dictionary<string, SX> { ["x"] = z, ["f"] = z * z },
        new Dictionary<string, object> { ["ipopt.linear_solver"] = "ma27", ["ipopt.print_level"] = 0L, ["print_time"] = false });
    s.call(new Dictionary<string, DM> { ["x0"] = 1 });
    return (bool)s.stats()["success"];
  }

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { Console.WriteLine("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }

    // Degree of interpolating polynomial
    int d = 3;

    // Choose collocation points
    double[] tau_root = ca.collocation_points(d, "radau").Prepend(0).ToArray();

    // Coefficients of the collocation equation
    var C = new double[d + 1, d + 1];

    // Coefficients of the continuity equation
    var D = new double[d + 1];

    // Coefficients of the quadrature function
    var F = new double[d + 1];

    // Construct polynomial basis
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      double[] p = { 1 };
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          double s = tau_root[j] - tau_root[r];
          p = polymul(p, [1 / s, -tau_root[r] / s]);
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D[j] = polyval(p, 1.0);

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      double[] pder = polyder(p);
      for (int r = 0; r < d + 1; r++) {
        C[j, r] = polyval(pder, tau_root[r]);
      }

      // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      double[] pint = polyint(p);
      F[j] = polyval(pint, 1.0);
    }

    // Control discretization
    int nk = 20;

    // End time
    double tf = 10.0;

    // Size of the finite elements
    double h = tf / nk;

    // All collocation time points
    var T = new double[nk, d + 1];
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < d + 1; j++)
        T[k, j] = h * (k + tau_root[j]);

    // Declare variables (use scalar graph)
    SX t = SX.sym("t");    // time
    SX u = SX.sym("u");    // control
    SX x = SX.sym("x", 2); // state

    // ODE rhs function and quadratures
    SX xdot = ca.vertcat((1 - x[1] * x[1]) * x[0] - x[1] + u, x[0]);
    SX qdot = x[0] * x[0] + x[1] * x[1] + u * u;
    var f = ca.Function("f", [t, x, u], [xdot, qdot]);

    // Control bounds
    double u_min = -0.75;
    double u_max = 1.0;
    double u_init = 0.0;

    // State bounds and initial guess
    double[] x_min = { -ca.inf, -ca.inf };
    double[] x_max = {  ca.inf,  ca.inf };
    double[] xi_min = { 0.0, 1.0 };
    double[] xi_max = { 0.0, 1.0 };
    double[] xf_min = { 0.0, 0.0 };
    double[] xf_max = { 0.0, 0.0 };
    double[] x_init = { 0.0, 0.0 };

    // Dimensions
    int nx = 2;
    int nu = 1;

    // Total number of variables
    int NX = nk * (d + 1) * nx; // Collocated states
    int NU = nk * nu;           // Parametrized controls
    int NXF = nx;               // Final state
    int NV = NX + NU + NXF;

    // NLP variable vector
    MX V = MX.sym("V", NV);

    // All variables with bounds and initial guess
    var vars_lb = new double[NV];
    var vars_ub = new double[NV];
    var vars_init = new double[NV];
    int offset = 0;

    // Get collocated states and parametrized control
    var X = new MX[nk + 1, d + 1];
    var U = new MX[nk];
    for (int k = 0; k < nk; k++) {
      // Collocated states
      for (int j = 0; j < d + 1; j++) {
        // Get the expression for the state vector
        X[k, j] = V[offset..(offset + nx)];

        // Add the initial condition
        Array.Copy(x_init, 0, vars_init, offset, nx);

        // Add bounds
        if (k == 0 && j == 0) {
          Array.Copy(xi_min, 0, vars_lb, offset, nx);
          Array.Copy(xi_max, 0, vars_ub, offset, nx);
        } else {
          Array.Copy(x_min, 0, vars_lb, offset, nx);
          Array.Copy(x_max, 0, vars_ub, offset, nx);
        }
        offset += nx;
      }

      // Parametrized controls
      U[k] = V[offset..(offset + nu)];
      Array.Fill(vars_lb, u_min, offset, nu);
      Array.Fill(vars_ub, u_max, offset, nu);
      Array.Fill(vars_init, u_init, offset, nu);
      offset += nu;
    }

    // State at end time
    X[nk, 0] = V[offset..(offset + nx)];
    Array.Copy(xf_min, 0, vars_lb, offset, nx);
    Array.Copy(xf_max, 0, vars_ub, offset, nx);
    Array.Copy(x_init, 0, vars_init, offset, nx);
    offset += nx;

    // Constraint function for the NLP
    var g = new List<MX>();
    int ng = 0;

    // Objective function
    MX J = 0;

    // For all finite elements
    for (int k = 0; k < nk; k++) {

      // For all collocation points
      for (int j = 1; j < d + 1; j++) {

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = 0;
        for (int r = 0; r < d + 1; r++) {
          xp_jk += C[r, j] * X[k, r];
        }

        // Add collocation equations to the NLP
        List<MX> fq = f.call([T[k, j], X[k, j], U[k]]);
        MX fk = fq[0], qk = fq[1];
        g.Add(h * fk - xp_jk);
        ng += nx; // equality constraints

        // Add contribution to objective
        J += F[j] * qk * h;
      }

      // Get an expression for the state at the end of the finite element
      MX xf_k = 0;
      for (int r = 0; r < d + 1; r++) {
        xf_k += D[r] * X[k, r];
      }

      // Add continuity equation to NLP
      g.Add(X[k + 1, 0] - xf_k);
      ng += nx;
    }

    // Concatenate constraints
    MX gc = ca.vertcat(g);

    // NLP
    var nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = J, ["g"] = gc };

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
    arg["lbg"] = new double[ng];
    arg["ubg"] = new double[ng];

    // Solve the problem
    var res = solver.call(arg);

    // Print the optimal cost
    Console.WriteLine($"optimal cost:  {(double)res["f"]}");

    // Retrieve the solution
    List<double> v_opt = res["x"].nonzeros();

    // Get values at the beginning of each finite element
    int stride = (d + 1) * nx + nu;
    double[] x0_opt = new double[nk + 1], x1_opt = new double[nk + 1], u_opt = new double[nk];
    for (int k = 0; k < nk + 1; k++) {
      x0_opt[k] = v_opt[k * stride];
      x1_opt[k] = v_opt[k * stride + 1];
    }
    for (int k = 0; k < nk; k++) u_opt[k] = v_opt[(d + 1) * nx + k * stride];
    DM tgrid = ca.linspace(0, tf, nk + 1);
    DM tgrid_u = ca.linspace(0, tf, nk);

    // plot: x0_opt, x1_opt over tgrid and u_opt (step) over tgrid_u, "Van der Pol optimization"
  }
}
