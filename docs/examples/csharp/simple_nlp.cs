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

public static class simple_nlp {
  public static void Main() {
    // Pick an NLP solver
    string MySolver = "ipopt";
    //string MySolver = "worhp";
    //string MySolver = "sqpmethod";
    if (!ca.has_nlpsol(MySolver)) { Console.WriteLine($"SKIP: nlpsol {MySolver} not available"); return; }

    // Declare variables
    SX x = SX.sym("x", 2);

    // Form the NLP
    SX f = ca.pow(x[0], 2) + ca.pow(x[1], 2); // objective
    SX g = x[0] + x[1] - 10;            // constraint
    var nlp = new Dictionary<string, SX> { ["x"] = x, ["f"] = f, ["g"] = g };

    // Solver options
    var opts = new Dictionary<string, object>();
    if (MySolver == "sqpmethod") {
      opts["qpsol"] = "qpoases";
      opts["qpsol_options"] = new Dictionary<string, object> { ["printLevel"] = "none" };
    }

    // Allocate a solver
    Function solver = ca.nlpsol("solver", MySolver, nlp, opts);

    // Solve the NLP
    var sol = solver.call(new Dictionary<string, DM> { ["lbg"] = 0 });

    // Print solution
    Console.WriteLine("-----");
    Console.WriteLine($"objective at solution =  {sol["f"]}");
    Console.WriteLine($"primal solution =  {sol["x"]}");
    Console.WriteLine($"dual solution (x) =  {sol["lam_x"]}");
    Console.WriteLine($"dual solution (g) =  {sol["lam_g"]}");
  }
}
