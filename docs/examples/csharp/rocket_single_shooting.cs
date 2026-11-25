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

public static class rocket_single_shooting {
  static bool sundials_integrator = true;
  static bool explicit_integrator = false;
  static bool lifted_newton = false;

  public static void Main() {
    if (!ca.has_integrator("cvodes") || !ca.has_integrator("idas") || !ca.has_integrator("rk")
        || !ca.has_nlpsol("ipopt") || !ca.has_nlpsol("scpgen")) {
      Console.WriteLine("SKIP: integrator cvodes/idas/rk or nlpsol ipopt/scpgen not available");
      return;
    }

    // Load integrators manually
    ca.load_integrator("cvodes");
    ca.load_integrator("idas");
    ca.load_integrator("rk");
    ca.load_nlpsol("ipopt");
    ca.load_nlpsol("scpgen");

    // Time length
    double T = 10.0;

    // Shooting length
    int nu = 20; // Number of control segments

    // Time horizon for integrator
    double t0 = 0;
    double tf = T / nu;

    // Initial position
    var X0 = new double[3];
    X0[0] = 0; // initial position
    X0[1] = 0; // initial speed
    X0[2] = 1; // initial mass

    // Time
    SX t = SX.sym("t");

    // Differential states
    SX s = SX.sym("s"), v = SX.sym("v"), m = SX.sym("m");
    SX x = ca.vertcat(s, v, m);

    // Control
    SX u = SX.sym("u");

    SX alpha = 0.05; // friction
    SX beta = 0.1; // fuel consumption rate

    // Differential equation
    SX rhs = ca.vertcat(v, (u - alpha * v * v) / m, -beta * u * u);

    // Initial conditions
    double[] x0 = { 0, 0, 1 };

    // DAE
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = u, ["t"] = t, ["ode"] = rhs };

    // Integrator options
    string plugin;
    var opts = new Dictionary<string, object>();
    if (sundials_integrator) {
      if (explicit_integrator) {
        // Explicit integrator (CVODES)
        plugin = "cvodes";
        // opts["exact_jacobian"] = true;
        // opts["linear_multistep_method"] = "bdf"; // adams or bdf
        // opts["nonlinear_solver_iteration"] = "newton"; // newton or functional
      } else {
        // Implicit integrator (IDAS)
        plugin = "idas";
        opts["calc_ic"] = false;
      }
      opts["fsens_err_con"] = true;
      opts["quad_err_con"] = true;
      opts["abstol"] = 1e-6;
      opts["reltol"] = 1e-6;
      //  opts["fsens_all_at_once"] = false;
      opts["steps_per_checkpoint"] = 100L; // BUG: Too low number causes segfaults
    } else {
      // An explicit Euler integrator
      plugin = "rk";
      opts["expand_f"] = true;
      opts["interpolation_order"] = 1L;
      opts["number_of_finite_elements"] = 1000L;
    }

    // Create integrator
    Function F = ca.integrator("integrator", plugin, dae, t0, tf, opts);

    // control for all segments
    MX U = MX.sym("U", nu);

    // Integrate over all intervals
    MX X = X0;
    for (int k = 0; k < nu; ++k) {
      // Integrate
      X = F.call(new Dictionary<string, MX> { ["x0"] = X, ["p"] = U[k] })["xf"];

      // Lift X
      if (lifted_newton) {
        // BINDING GAP: MX lift(x, x_guess) is not wrapped in the C# binding (as in Java); it should read
        // X = lift(X, X);
        Console.WriteLine("SKIP: lift(MX, MX) is not wrapped in the C# binding (as in Java)");
        return;
      }
    }

    // Objective function
    MX J = ca.dot(U, U);

    // Terminal constraints
    MX G = ca.vertcat(X[0], X[1]);

    // Create the NLP
    var nlp = new Dictionary<string, MX> { ["x"] = U, ["f"] = J, ["g"] = G };

    // NLP solver options
    var solver_opts = new Dictionary<string, object>();
    string solver_name;
    if (lifted_newton) {
      solver_name = "scpgen";
      solver_opts["verbose"] = true;
      solver_opts["regularize"] = false;
      solver_opts["max_iter_ls"] = 1L;
      solver_opts["max_iter"] = 100L;
      solver_opts["qpsol"] = "nlp"; // Use IPOPT as QP solver
      var ipopt_options = new Dictionary<string, object> {
        ["tol"] = 1e-12,
        ["print_level"] = 0L,
        ["print_time"] = false };
      solver_opts["qpsol_options"] = new Dictionary<string, object> {
        ["nlpsol_options"] = new Dictionary<string, object> { ["nlpsol"] = "ipopt" } };
    } else {
      solver_name = "ipopt";
      solver_opts["ipopt.tol"] = 1e-10;
      solver_opts["ipopt.hessian_approximation"] = "limited-memory";
    }

    // Create NLP solver
    Function solver = ca.nlpsol("nlpsol", solver_name, nlp, solver_opts);

    // Bounds on u and initial condition
    double[] umin = new double[nu], umax = new double[nu], u0 = new double[nu];
    Array.Fill(umin, -10);
    Array.Fill(umax, 10);
    Array.Fill(u0, 0.4);

    // Bounds on g
    double[] gmin = { 10, 0 }, gmax = { 10, 0 };

    // Solve NLP
    DM Usol = solver.call(new Dictionary<string, DM> {
      ["lbx"] = umin, ["ubx"] = umax, ["x0"] = u0,
      ["lbg"] = gmin, ["ubg"] = gmax })["x"];
    Console.WriteLine($"optimal solution: {Usol}");
  }
}
