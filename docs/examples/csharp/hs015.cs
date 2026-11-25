//
//     MIT No Attribution
//
//     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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

/*
Solve the hs015 problem, formulated as the NLP:

minimize     100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
subject to   x[0]*x[1] >= 1
             x[0] + x[1]**2 >= 0
             x[0] <= 1/2

David Kiessling, 2026
*/

public static class hs015 {
  public static void Main() {
    if (!ca.has_nlpsol("uno")) { Console.WriteLine("SKIP: nlpsol uno not available"); return; }

    // Declare variable x
    MX x = MX.sym("x", 2);

    // Formulate the NLP
    MX f = 100 * ca.pow(x[1] - ca.pow(x[0], 2), 2) + ca.pow(1 - x[0], 2);
    MX g = MX.zeros(2, 1);
    g[0] = x[0] * x[1];
    g[1] = x[0] + ca.pow(x[1], 2);

    DM ubx = new double[] { 1.0 / 2, ca.inf };
    DM lbx = new double[] { -ca.inf, -ca.inf };
    DM lbg = new double[] { 1, 0 };
    DM ubg = new double[] { ca.inf, ca.inf };
    var nlp = new Dictionary<string, MX> { ["x"] = x, ["f"] = f, ["g"] = g };

    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "uno", nlp, new Dictionary<string, object> {
      ["uno"] = new Dictionary<string, object> { ["preset"] = "ipopt" } });

    // Solve the problem
    var res = solver.call(new Dictionary<string, DM> {
      ["x0"] = new double[] { -2.0, 1.0 },
      ["ubg"] = ubg,
      ["lbg"] = lbg,
      ["ubx"] = ubx });

    // Print solution
    Console.WriteLine();
    Console.WriteLine($"{"Optimal cost:",50}  {res["f"]}");
    Console.WriteLine($"{"Primal solution:",50}  {res["x"]}");
    Console.WriteLine($"{"Dual solution (simple bounds):",50}  {res["lam_x"]}");
    Console.WriteLine($"{"Dual solution (nonlinear bounds):",50}  {res["lam_g"]}");
  }
}
