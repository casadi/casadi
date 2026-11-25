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
Solve a basic mpcc with ccopt

minimize     (x-2)^2 + (y-2)^2
subject to   0 <= x _|_ y >= 0

*/

public static class basic_mpcc {
  public static void Main() {
    if (!ca.has_nlpsol("ccopt")) { Console.WriteLine("SKIP: nlpsol ccopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");

    // Formulate the NLP
    SX f = ca.pow(x - 2, 2) + ca.pow(y - 2, 2);
    var nlp = new Dictionary<string, SX> { ["x"] = ca.vertcat(x, y),
                                           ["f"] = f };

    var solver_opts = new Dictionary<string, object> {
      //["madnlp.print_level"] = 1L,
      ["cc_pairs"] = new List<List<long>> { new List<long> { 0, 1 } },
      ["madnlp.bound_relax_factor"] = 0.0 };
    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "ccopt", nlp, solver_opts);

    // Solve the Rosenbrock problem
    var res = solver.call(new Dictionary<string, DM> {
      ["x0"] = new double[] { 2.5, 3.0 },
      ["lbx"] = new double[] { 0.0, 0.0 } });

    //  Print solution
    Console.WriteLine($"Optimal cost:                     {(double)res["f"]}");
    Console.WriteLine($"Primal solution:                  {res["x"]}");
    Console.WriteLine($"Dual solution (simple bounds):    {res["lam_x"]}");
  }
}
