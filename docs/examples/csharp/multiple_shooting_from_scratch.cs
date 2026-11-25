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

/** Writing a multiple shooting code from scratch

  This example demonstrates how to write a simple, yet powerful multiple shooting code from
  scratch using CasADi. For clarity, the code below uses a simple formulation with only states
  and controls, and no path constrants. It relies on CasADi machinery to keep track of sparsity,
  formulate ODE sensitivity equations and build up the Jacobian of the NLP constraint function.

  By extending the code below, it should be possible for a user to solve ever more complex
  problems. For example, one can easily make a multi-stage formulation by simply allocating
  another integrator instance and use the two integrators instances for different shooting nodes.
  By replacing the explicit CVodes integrator with the fully-implicit IDAS integrator, one can
  solve optimal control problems in differential-algebraic equations.

  Joel Andersson, 2011-2012
*/

public static class multiple_shooting_from_scratch {
  public static void Main() {
    if (!ca.has_integrator("cvodes")) { Console.WriteLine("SKIP: integrator cvodes not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Declare variables
    SX u = SX.sym("u"); // control
    SX r = SX.sym("r"), s = SX.sym("s"); // states
    SX x = ca.vertcat(r, s);

    // Number of differential states
    long nx = x.size1();

    // Number of controls
    long nu = u.size1();

    // Bounds and initial guess for the control
    double[] u_min  = { -0.75 };
    double[] u_max  = {  1.0  };
    double[] u_init = {  0.0  };

    // Bounds and initial guess for the state
    double[] x0_min = {   0.0,  1.0 };
    double[] x0_max = {   0.0,  1.0 };
    double[] x_min  = {  -ca.inf, -ca.inf };
    double[] x_max  = {   ca.inf,  ca.inf };
    double[] xf_min = {   0.0,  0.0 };
    double[] xf_max = {   0.0,  0.0 };
    double[] x_init = {   0.0,  0.0 };

    // Final time
    double tf = 20.0;

    // Number of shooting nodes
    int ns = 50;

    // ODE right hand side and quadrature
    SX ode = ca.vertcat((1 - s * s) * r - s + u, r);
    SX quad = r * r + s * s + u * u;
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = u, ["ode"] = ode, ["quad"] = quad };

    // Create an integrator (CVodes)
    Function F = ca.integrator("integrator", "cvodes", dae, 0, tf / ns);

    // Total number of NLP variables
    long NV = nx * (ns + 1) + nu * ns;

    // Declare variable vector for the NLP
    MX V = MX.sym("V", NV);

    // NLP variable bounds and initial guess
    List<double> v_min = new(), v_max = new(), v_init = new();

    // Offset in V
    long offset = 0;

    // State at each shooting node and control for each shooting interval
    List<MX> X = new(), U = new();
    for (int k = 0; k < ns; ++k) {
      // Local state
      X.Add(V.get_nz(ca.Slice(offset, offset + nx)));
      if (k == 0) {
        v_min.AddRange(x0_min);
        v_max.AddRange(x0_max);
      } else {
        v_min.AddRange(x_min);
        v_max.AddRange(x_max);
      }
      v_init.AddRange(x_init);
      offset += nx;

      // Local control
      U.Add(V.get_nz(ca.Slice(offset, offset + nu)));
      v_min.AddRange(u_min);
      v_max.AddRange(u_max);
      v_init.AddRange(u_init);
      offset += nu;
    }

    // State at end
    X.Add(V.get_nz(ca.Slice(offset, offset + nx)));
    v_min.AddRange(xf_min);
    v_max.AddRange(xf_max);
    v_init.AddRange(x_init);
    offset += nx;

    // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
    if (offset != NV) throw new Exception("");

    // Objective function
    MX J = 0;

    //Constraint function and bounds
    var g = new List<MX>();

    // Loop over shooting nodes
    for (int k = 0; k < ns; ++k) {
      // Create an evaluation node
      var I_out = F.call(new Dictionary<string, MX> { ["x0"] = X[k], ["p"] = U[k] });

      // Save continuity constraints
      g.Add(I_out["xf"] - X[k + 1]);

      // Add objective function contribution
      J += I_out["qf"];
    }

    // NLP
    var nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = J, ["g"] = ca.vertcat(g) };

    // Set options
    var opts = new Dictionary<string, object>();
    opts["ipopt.tol"] = 1e-5;
    opts["ipopt.max_iter"] = 100L;

    // Create an NLP solver and buffers
    Function solver = ca.nlpsol("nlpsol", "ipopt", nlp, opts);
    var arg = new Dictionary<string, DM>();

    // Bounds and initial guess
    arg["lbx"] = v_min.ToArray();
    arg["ubx"] = v_max.ToArray();
    arg["lbg"] = 0;
    arg["ubg"] = 0;
    arg["x0"] = v_init.ToArray();

    // Solve the problem
    var res = solver.call(arg);

    // Optimal solution of the NLP
    List<double> V_opt = res["x"].nonzeros();

    // Get the optimal state trajectory
    double[] r_opt = new double[ns + 1], s_opt = new double[ns + 1];
    for (int i = 0; i <= ns; ++i) {
      r_opt[i] = V_opt[(int)(i * (nx + 1))];
      s_opt[i] = V_opt[(int)(1 + i * (nx + 1))];
    }
    Console.WriteLine("r_opt = ");
    Console.WriteLine((DM)r_opt);
    Console.WriteLine("s_opt = ");
    Console.WriteLine((DM)s_opt);

    // Get the optimal control
    var u_opt = new double[ns];
    for (int i = 0; i < ns; ++i) {
      u_opt[i] = V_opt[(int)(nx + i * (nx + 1))];
    }
    Console.WriteLine("u_opt = ");
    Console.WriteLine((DM)u_opt);
  }
}
