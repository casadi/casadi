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

/** Solves the following NLP with block structured Hessian:
    minimize      x^2 - 0.5*y^2
    subject to    x - y = 0

    Joel Andersson, 2016
*/
public static class test_blocksqp {
  public static void Main() {
    if (!ca.has_nlpsol("blocksqp")) { Console.WriteLine("SKIP: nlpsol blocksqp not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");

    // Formulate the NLP
    SX f = ca.pow(x, 2) - 0.5 * ca.pow(y, 2);
    SX g = x - y;
    var nlp = new Dictionary<string, SX> { ["x"] = ca.vertcat(x, y), ["f"] = f, ["g"] = g };

    // Create an NLP solver
    var opts = new Dictionary<string, object> {
      ["opttol"] = 1.0e-12,
      ["nlinfeastol"] = 1.0e-12,
      ["globalization"] = 0L,
      ["hess_update"] = 0L,
      ["hess_scaling"] = 0L,
      ["fallback_scaling"] = 0L,
      ["hess_lim_mem"] = 0L,
      ["max_consec_skipped_updates"] = 200L,
      ["block_hess"] = 0L,
      ["which_second_derv"] = 0L,
      ["schur"] = false,
    };
    Function solver = ca.nlpsol("solver", "blocksqp", nlp, opts);

    // Solve the Rosenbrock problem
    var res = solver.call(new Dictionary<string, DM> {
      ["x0"] = new double[] { 10, 10 },
      ["lbg"] = 0,
      ["ubg"] = 0 });

    //  Print solution
    Console.WriteLine($"Optimal cost:                     {(double)res["f"]}");
    Console.WriteLine($"Primal solution:                  {res["x"]}");
    Console.WriteLine($"Dual solution (simple bounds):    {res["lam_x"]}");
    Console.WriteLine($"Dual solution (nonlinear bounds): {res["lam_g"]}");
  }
}
