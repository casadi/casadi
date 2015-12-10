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

  /* Initialize matrix A. */
  int ncol = 5, nrow = 5;
  int nnz = 12;
  
  vector<int> colind(ncol+1);
  vector<int> row(nnz);
  
  // Sparsity pattern
  row[0] = 0; row[1] = 1; row[2] = 4; row[3] = 1;
  row[4] = 2; row[5] = 4; row[6] = 0; row[7] = 2;
  row[8] = 0; row[9] = 3; row[10]= 3; row[11]= 4;
  colind[0] = 0; colind[1] = 3; colind[2] = 6; colind[3] = 8; colind[4] = 10; colind[5] = 12;
  Sparsity spA(nrow,ncol,colind,row);
  
  // Create a solver instance
  Function linear_solver = linsol("linear_solver", "csparse", spA, 1);
    
  // Pass Non-zero elements
  double   s, u, p, e, r, l;
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;

  vector<double> a(nnz);
  a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
  a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
  
  // Right hand side
  vector<double> b(ncol,1.0);

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
