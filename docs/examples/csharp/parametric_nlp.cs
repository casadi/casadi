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

/**
 *  Example program demonstrating parametric NLPs in CasADi
 *  Note that there is currently no support for parametric sensitivities via this feature (although it would make a lot of sense).
 *  For parametric sensitivities, see the parametric_sensitivities.cpp example which calculates sensitivitied via the sIPOPT extension
 *  to IPOPT.
 *
 *  Joel Andersson, KU Leuven 2012
 */

public static class parametric_nlp {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    /* Test problem (Ganesh & Biegler, A reduced Hessian strategy for sensitivity analysis of optimal flowsheets, AIChE 33, 1987, pp. 282-296)
     *
     *    min     x1^2 + x2^2 + x3^2
     *    s.t.    6*x1 + 3&x2 + 2*x3 - pi = 0
     *            p2*x1 + x2 - x3 - 1 = 0
     *            x1, x2, x3 >= 0
     *
     */

    // Optimization variables
    SX x = SX.sym("x", 3);

    // Parameters
    SX p = SX.sym("p", 2);

    // Objective
    SX f = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

    // Constraints
    SX g = ca.vertcat(
        6 * x[0] + 3 * x[1] + 2 * x[2] - p[0],
        p[1] * x[0] + x[1] - x[2] - 1
    );

    // Initial guess and bounds for the optimization variables
    double[] x0  = { 0.15, 0.15, 0.00 };
    double[] lbx = { 0.00, 0.00, 0.00 };
    double[] ubx = { ca.inf, ca.inf, ca.inf };

    // Nonlinear bounds
    double[] lbg = { 0.00, 0.00 };
    double[] ubg = { 0.00, 0.00 };

    // Original parameter values
    double[] p0  = { 5.00, 1.00 };

    // NLP
    var nlp = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["f"] = f, ["g"] = g };

    // Create NLP solver and buffers
    Function solver = ca.nlpsol("solver", "ipopt", nlp);
    var arg = new Dictionary<string, DM>();

    // Solve the NLP
    arg["lbx"] = lbx;
    arg["ubx"] = ubx;
    arg["lbg"] = lbg;
    arg["ubg"] = ubg;
    arg["x0"] = x0;
    arg["p"] = p0;
    var res = solver.call(arg);

    // Print the solution
    Console.WriteLine("-----");
    Console.WriteLine($"Optimal solution for p = {arg["p"]}:");
    Console.WriteLine($"{"Objective: ",30}{res["f"]}");
    Console.WriteLine($"{"Primal solution: ",30}{res["x"]}");
    Console.WriteLine($"{"Dual solution (x): ",30}{res["lam_x"]}");
    Console.WriteLine($"{"Dual solution (g): ",30}{res["lam_g"]}");

    // Change the parameter and resolve
    arg["p"] = 4.5;
    res = solver.call(arg);

    // Print the new solution
    Console.WriteLine("-----");
    Console.WriteLine($"Optimal solution for p = {arg["p"]}:");
    Console.WriteLine($"{"Objective: ",30}{res["f"]}");
    Console.WriteLine($"{"Primal solution: ",30}{res["x"]}");
    Console.WriteLine($"{"Dual solution (x): ",30}{res["lam_x"]}");
    Console.WriteLine($"{"Dual solution (g): ",30}{res["lam_g"]}");
  }
}
