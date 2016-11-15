/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
File superlu.c from the CSparse example collection
Joel Andersson, K.U. Leuven, 2010
*/

#include "casadi/casadi.hpp"

using namespace casadi;
using namespace std;

int main(int argc, char *argv[])
{

  // A
  int ncol = 5, nrow = 5;
  vector<int> colind = {0, 3, 6, 8, 10, 12};
  vector<int> row = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4};
  vector<double> nz = {19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};
  DM A(Sparsity(nrow, ncol, colind, row), nz);

  // Right hand side
  DM b = DM::ones(ncol);

  // Type of linear systems
  enum SymType {UNSYM, SYM, PD};

  // All Linear solvers to be tested
  struct Test {
    string solver;
    SymType type;
  };
  vector<Test> tests;
  tests.push_back({"csparse", UNSYM});
  tests.push_back({"csparsecholesky", PD});
  tests.push_back({"lapacklu", UNSYM});
  tests.push_back({"lapackqr", UNSYM});
  tests.push_back({"ma27", SYM});
  tests.push_back({"symbolicqr", UNSYM});

  // Test all combinations
  for (auto s : {UNSYM, SYM, PD}) {
    DM A_test;
    switch (s) {
    case UNSYM:
      cout << "Unsymmetric linear system" << endl;
      A_test = A;
      break;
    case SYM:
      cout << "Symmetric linear system" << endl;
      A_test = A + A.T();
      break;
    case PD:
      cout << "Positive definite linear system" << endl;
      A_test = mtimes(A.T(), A);
      break;
    }
    for (auto t : tests) {
      if (t.type > s) continue; // Cannot be solved
      if (!Linsol::has_plugin(t.solver)) {
        cout << t.solver << " not available" << endl;
        continue;
      }

      // Create a solver instance
      Linsol F("F", t.solver);

      // Solve
      F.reset(A_test.sparsity());
      F.pivoting(A_test.ptr());
      F.factorize(A_test.ptr());
      DM x = densify(b);
      F.solve(x.ptr(), x.size2());

      // Print the solution
      cout << "solution: " << x << " (" <<  t.solver << ")" << endl;
    }
  }

  return 0;
}
