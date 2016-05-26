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

  // Sparsity pattern of A
  int ncol = 5, nrow = 5;
  vector<int> colind = {0, 3, 6, 8, 10, 12};
  vector<int> row = {0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4};
  Sparsity spA(nrow, ncol, colind, row);

  // Create a solver instance
  Function linear_solver = linsol("linear_solver", "csparse", spA, 1);

  // Nonzeros of A
  vector<double> a = {19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18};

  // Right hand side
  vector<double> b(ncol, 1.0);

  // Solve
  vector<double> x;
  linear_solver({{"A", a}, {"B", b}}, {{"X", &x}});

  // Print the solution
  cout << "solution            = " << x << endl;

  // Embed in an MX graph
  MX A = MX::sym("A", spA);
  MX B = MX::sym("B", ncol);
  MX X = linear_solver.linsol_solve(A, B);
  Function F("F", {A, B}, {X}, {"A", "B"}, {"X"});

  // Solve
  F({{"A", a}, {"B", b}}, {{"X", &x}});

  // Print the solution
  cout << "solution (embedded) = " << x << endl;

  return 0;
}
