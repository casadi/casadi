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
import org.casadi.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
File superlu.c from the CSparse example collection
Joel Andersson, KU Leuven, 2010
*/

public class test_linsol {
  // Type of linear systems
  enum SymType {UNSYM, SYM, PD}

  // All Linear solvers to be tested
  record Test(String solver, SymType type) {}

  public static void main(String[] args) {
    // A
    int ncol = 5, nrow = 5;
    long[] colind = {0, 3, 6, 8, 10, 12};
    long[] row = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4};
    double[] nz = {19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};
    DM A = new DM(new Sparsity(nrow, ncol, colind, row), nz);

    // Right hand side
    DM b = DM.ones(ncol);

    List<Test> tests = new ArrayList<>();
    tests.add(new Test("csparse", SymType.UNSYM));
    tests.add(new Test("csparsecholesky", SymType.PD));
    tests.add(new Test("lapacklu", SymType.UNSYM));
    tests.add(new Test("lapackqr", SymType.UNSYM));
    tests.add(new Test("ma27", SymType.SYM));
    tests.add(new Test("mumps", SymType.UNSYM));
    tests.add(new Test("symbolicqr", SymType.UNSYM));
    tests.add(new Test("qr", SymType.UNSYM));
    tests.add(new Test("ldl", SymType.SYM));
    tests.add(new Test("lsqr", SymType.UNSYM));

    // Test all combinations
    for (SymType s : SymType.values()) {
      DM A_test = null;
      switch (s) {
      case UNSYM:
        System.out.println("Unsymmetric linear system");
        A_test = A;
        break;
      case SYM:
        System.out.println("Symmetric linear system");
        A_test = A.plus(A.T());
        break;
      case PD:
        System.out.println("Positive definite linear system");
        A_test = ca.mtimes(A.T(), A);
        break;
      }
      for (Test t : tests) {
        if (t.type().compareTo(s) > 0) continue; // Cannot be solved
        if (!Linsol.has_plugin(t.solver())) {
          System.out.println(t.solver() + " not available");
          continue;
        }

        // Solver specific options
        Map<String, Object> opts = new HashMap<>();
        if (t.solver().equals("mumps")) {
          opts.put("symmetric", s == SymType.SYM || s == SymType.PD);
          opts.put("posdef", s == SymType.PD);
        }

        // Create a solver instance
        Linsol F = new Linsol("F", t.solver(), A_test.sparsity(), opts);

        // Solve (failures are thrown as CasadiException)
        F.sfact(A_test);
        F.nfact(A_test);
        DM x = F.solve(A_test, ca.densify(b));

        // Print the solution
        System.out.println("solution: " + x + " (" + t.solver() + ")");
      }
    }
  }
}
