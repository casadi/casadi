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

/*
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2 + 100*z^2
subject to   z+(1-x)^2-y == 0

Joel Andersson, 2015
*/

public static class rosenbrock {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");

    // Formulate the NLP
    SX f = ca.pow(x, 2) + 100 * ca.pow(z, 2);
    SX g = z + ca.pow(1 - x, 2) - y;
    var nlp = new Dictionary<string, SX> { ["x"] = ca.vertcat(x, y, z), ["f"] = f, ["g"] = g };

    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp);

    // Solve the Rosenbrock problem
    var res = solver.call(new Dictionary<string, DM> {
      ["x0"] = new double[] { 2.5, 3.0, 0.75 },
      ["ubg"] = 0,
      ["lbg"] = 0 });

    // Print solution
    Console.WriteLine();
    Console.WriteLine($"{"Optimal cost:",50}  {res["f"]}");
    Console.WriteLine($"{"Primal solution:",50}  {res["x"]}");
    Console.WriteLine($"{"Dual solution (simple bounds):",50}  {res["lam_x"]}");
    Console.WriteLine($"{"Dual solution (nonlinear bounds):",50}  {res["lam_g"]}");
  }
}
