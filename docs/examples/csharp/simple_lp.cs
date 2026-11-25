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

// minimize    3x + 4y
// subject to  x + 2y <= 14
//            3x -  y >= 0
//             x -  y <= 2

public static class simple_lp {
  public static void Main() {
    if (!ca.has_conic("qpoases")) { Console.WriteLine("SKIP: conic qpoases not available"); return; }

    // Sparsity of the LP linear term
    Sparsity A = Sparsity.dense(3, 2);

    // Create solver
    Function solver = ca.conic("solver", "qpoases", new Dictionary<string, Sparsity> { ["a"] = A });
    //Function solver = conic("solver", "clp", new Dictionary<string, Sparsity> { ["a"] = A }); // Use clp

    DM g = new double[] { 3, 4 };
    DM a = new double[,] { { 1, 2 }, { 3, -1 }, { 1, -1 } };
    DM lba = new double[] { -ca.inf, 0, -ca.inf };
    DM uba = new double[] { 14, ca.inf, 2 };

    var sol = solver.call(new Dictionary<string, DM> { ["g"] = g, ["a"] = a, ["lba"] = lba, ["uba"] = uba });
    Console.WriteLine("{" + string.Join(", ", sol.Select(kv => $"{kv.Key}={kv.Value}")) + "}");
  }
}
