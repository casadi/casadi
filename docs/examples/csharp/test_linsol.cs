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
File superlu.c from the CSparse example collection
Joel Andersson, KU Leuven, 2010
*/

public static class test_linsol {
  // Type of linear systems
  enum SymType { UNSYM, SYM, PD }

  // All Linear solvers to be tested
  record Test(string solver, SymType type);

  public static void Main() {
    // A
    int ncol = 5, nrow = 5;
    long[] colind = { 0, 3, 6, 8, 10, 12 };
    long[] row = { 0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4 };
    double[] nz = { 19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18 };
    var A = ca.DM(ca.Sparsity(nrow, ncol, colind, row), nz);

    // Right hand side
    DM b = DM.ones(ncol);

    var tests = new List<Test> {
      new("csparse", SymType.UNSYM),
      new("csparsecholesky", SymType.PD),
      new("lapacklu", SymType.UNSYM),
      new("lapackqr", SymType.UNSYM),
      new("ma27", SymType.SYM),
      new("mumps", SymType.UNSYM),
      new("symbolicqr", SymType.UNSYM),
      new("qr", SymType.UNSYM),
      new("ldl", SymType.SYM),
      new("lsqr", SymType.UNSYM),
    };

    // Test all combinations
    foreach (SymType s in Enum.GetValues<SymType>()) {
      DM A_test = null;
      switch (s) {
      case SymType.UNSYM:
        Console.WriteLine("Unsymmetric linear system");
        A_test = A;
        break;
      case SymType.SYM:
        Console.WriteLine("Symmetric linear system");
        A_test = A + A.T();
        break;
      case SymType.PD:
        Console.WriteLine("Positive definite linear system");
        A_test = ca.mtimes(A.T(), A);
        break;
      }
      foreach (Test t in tests) {
        if (t.type > s) continue; // Cannot be solved
        if (!Linsol.has_plugin(t.solver)) {
          Console.WriteLine($"{t.solver} not available");
          continue;
        }

        // Solver specific options
        var opts = new Dictionary<string, object>();
        if (t.solver == "mumps") {
          opts["symmetric"] = s == SymType.SYM || s == SymType.PD;
          opts["posdef"] = s == SymType.PD;
        }

        // Create a solver instance
        var F = ca.Linsol("F", t.solver, A_test.sparsity(), opts);

        // Solve (failures are thrown as SwigException)
        F.sfact(A_test);
        F.nfact(A_test);
        DM x = F.solve(A_test, ca.densify(b));

        // Print the solution
        Console.WriteLine($"solution: {x} ({t.solver})");
      }
    }
  }
}
