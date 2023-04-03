/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

/**
File superlu.c from the CSparse example collection
Joel Andersson, KU Leuven, 2010
*/

#include "casadi/casadi.hpp"

using namespace casadi;

int main(int argc, char *argv[])
{

  // A
  int ncol = 5, nrow = 5;
  std::vector<casadi_int> colind = {0, 3, 6, 8, 10, 12};
  std::vector<casadi_int> row = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4};
  std::vector<double> nz = {19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};
  DM A(Sparsity(nrow, ncol, colind, row), nz);

  // Right hand side
  DM b = DM::ones(ncol);

  // Type of linear systems
  enum SymType {UNSYM, SYM, PD};

  // All Linear solvers to be tested
  struct Test {
    std::string solver;
    SymType type;
  };
  std::vector<Test> tests;
  tests.push_back({"csparse", UNSYM});
  tests.push_back({"csparsecholesky", PD});
  tests.push_back({"lapacklu", UNSYM});
  tests.push_back({"lapackqr", UNSYM});
  tests.push_back({"ma27", SYM});
  tests.push_back({"mumps", UNSYM});
  tests.push_back({"symbolicqr", UNSYM});
  tests.push_back({"qr", UNSYM});
  tests.push_back({"ldl", SYM});
  tests.push_back({"lsqr", UNSYM});

  // Test all combinations
  for (auto s : {UNSYM, SYM, PD}) {
    DM A_test;
    switch (s) {
    case UNSYM:
      std::cout << "Unsymmetric linear system" << std::endl;
      A_test = A;
      break;
    case SYM:
      std::cout << "Symmetric linear system" << std::endl;
      A_test = A + A.T();
      break;
    case PD:
      std::cout << "Positive definite linear system" << std::endl;
      A_test = mtimes(A.T(), A);
      break;
    }
    for (auto t : tests) {
      if (t.type > s) continue; // Cannot be solved
      if (!Linsol::has_plugin(t.solver)) {
        std::cout << t.solver << " not available" << std::endl;
        continue;
      }

      // Solver specific options
      Dict opts;
      if (t.solver == "mumps") {
        opts["symmetric"] = s == SYM || s == PD;
        opts["posdef"] = s == PD;
      }

      // Create a solver instance
      Linsol F("F", t.solver, A_test.sparsity(), opts);

      // Solve
      if (F.sfact(A_test.ptr())) casadi_error("'sfact' failed");
      if (F.nfact(A_test.ptr())) casadi_error("'nfact' failed");
      DM x = densify(b);
      if (F.solve(A_test.ptr(), x.ptr(), x.size2())) casadi_error("'solve' failed");

      // Print the solution
      std::cout << "solution: " << x << " (" <<  t.solver << ")" << std::endl;
    }
  }

  return 0;
}
