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

public static class vdp_indirect_multiple_shooting {
  public static void Main() {
    if (!ca.has_integrator("cvodes")) { Console.WriteLine("SKIP: integrator cvodes not available"); return; }
    if (!ca.has_rootfinder("nlpsol") || !ca.has_nlpsol("ipopt")) {
      Console.WriteLine("SKIP: rootfinder nlpsol with ipopt not available");
      return;
    }

    // Declare variables (use simple, efficient DAG)
    SX x0 = SX.sym("x0"); SX x1 = SX.sym("x1");
    SX x = ca.vertcat(x0, x1);

    // Control
    SX u = SX.sym("u");

    // ODE right hand side
    SX xdot = ca.vertcat((1 - x1 * x1) * x0 - x1 + u, x0);

    // Lagrangian function
    SX L = x0 * x0 + x1 * x1 + u * u;

    // Costate
    SX lam = SX.sym("lam", 2);

    // Hamiltonian function
    SX H = ca.dot(lam, xdot) + L;

    // Costate equations
    SX ldot = -ca.gradient(H, x);

    //# The control must minimize the Hamiltonian, which is:
    Console.WriteLine($"Hamiltonian:  {H}");

    // H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
    SX p = ca.gradient(H, u);     // this gives us 2*u + p
    p = ca.substitute(p, u, 0); // replace u with zero: gives us p

    // H's unconstrained minimizer is: u = -p/2
    SX u_opt = -p / 2;

    // We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
    u_opt = ca.fmin(u_opt, 1.0);
    u_opt = ca.fmax(u_opt, -0.75);
    Console.WriteLine($"optimal control:  {u_opt}");

    // Augment f with lam_dot and subtitute in the value for the optimal control
    SX f = ca.vertcat(xdot, ldot);
    f = ca.substitute(f, u, u_opt);

    // Function for obtaining the optimal control from the augmented state
    var u_fcn = ca.Function("ufcn", [ca.vertcat(x, lam)], [u_opt]);

    // Formulate the DAE
    var dae = new Dictionary<string, SX> { ["x"] = ca.vertcat(x, lam), ["ode"] = f };

    // Augmented DAE state dimension
    int nX = 4;

    // End time
    double tf = 10.0;

    // Number of shooting nodes
    int num_nodes = 20;

    // Create an integrator (CVodes)
    var iopts = new Dictionary<string, object>();
    iopts["abstol"] = 1e-8; // abs. tolerance
    iopts["reltol"] = 1e-8; // rel. tolerance
    Function I = ca.integrator("I", "cvodes", dae, 0, tf / num_nodes, iopts);

    // Variables for the states at each shooting node
    MX X = MX.sym("X", nX, num_nodes + 1);

    // Formulate the root finding problem
    var G = new List<MX>();
    G.Add(X[0..2, 0] - new double[] { 0, 1 }); // states fixed, costates free at initial time
    for (int k = 0; k < num_nodes; k++) {
      MX XF = I.call(new Dictionary<string, MX> { ["x0"] = X[.., k] })["xf"];
      G.Add(XF - X[.., k + 1]);
    }
    G.Add(X[2..nX, num_nodes] - new double[] { 0, 0 }); // costates fixed, states free at final time

    // Terminal constraints: lam = 0
    var rfp = new Dictionary<string, MX> { ["x"] = ca.vec(X), ["g"] = ca.vertcat(G) };

    // Select a solver for the root-finding problem
    string Solver = "nlpsol";
    //string Solver = "newton";
    //string Solver = "kinsol";

    // Solver options
    var opts = new Dictionary<string, object>();
    if (Solver == "nlpsol") {
      opts["nlpsol"] = "ipopt";
      opts["nlpsol_options"] = new Dictionary<string, object> { ["ipopt.hessian_approximation"] = "limited-memory" };
    } else if (Solver == "newton") {
      opts["linear_solver"] = "csparse";
    } else if (Solver == "kinsol") {
      opts["linear_solver_type"] = "user_defined";
      opts["linear_solver"] = "csparse";
      opts["max_iter"] = 1000L;
    }

    // Allocate a solver
    Function solver = ca.rootfinder("solver", Solver, rfp, opts);

    // Solve the problem
    DM X_sol = solver.call(new Dictionary<string, DM> { ["x0"] = 0 })["x"];
    Console.WriteLine(X_sol);

    // Time grid for visualization
    double[] tgrid = ca.linspace(0, tf, 100).nonzeros().ToArray();

    // Simulator to get optimal state and control trajectories
    Function simulator = ca.integrator("simulator", "cvodes", dae, 0, tgrid);

    // Simulate to get the trajectories
    DM sol = simulator.call(new Dictionary<string, DM> { ["x0"] = X_sol[0..4] })["xf"];

    // Calculate the optimal control
    DM u_opt_traj = u_fcn.Call(sol);

    // plot: sol[0,:], sol[1,:] and u_opt versus tgrid, 'Van der Pol optimization - indirect multiple shooting'

    Console.WriteLine(sol[0, ..]);
  }
}
