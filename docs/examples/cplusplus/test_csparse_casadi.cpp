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
  LinearSolver linear_solver("csparse", spA);
    
  // Initialize
  linear_solver.init();

  // Pass Non-zero elements
  double   s, u, p, e, r, l;
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;

  vector<double> val(nnz);
  val[0] = s; val[1] = l; val[2] = l; val[3] = u; val[4] = l; val[5] = l;
  val[6] = u; val[7] = p; val[8] = u; val[9] = e; val[10]= u; val[11]= r;
  
  // Right hand side
  vector<double> rhs(ncol,1.0);
  
  // Transpose?
  bool tr = false;

  // Solve
  linear_solver.setInput(val,"A");
  linear_solver.setInput(rhs,"B");
  linear_solver.prepare();
  linear_solver.solve(tr);
  
  // Print the solution
  cout << "solution = " << linear_solver.output("X") << endl;

  // Embed in an MX graph
  MX A = MX::sym("A",spA);
  MX B = MX::sym("B",ncol,1);
  MX X = linear_solver.solve(A,B,tr);
  MXFunction F(linsolIn("A",A,"B",B),linsolOut("X",X));
  F.init();

  // Solve
  F.setInput(val,"A");
  F.setInput(rhs,"B");
  F.evaluate();
  
  // Print the solution
  cout << "solution (embedded) = " << F.output("X") << endl;
  
  // Preturb the linear solver
  double t = 0.01;
  DMatrix x_unpreturbed = F.output("X");
  F.input("A")(3,2)   += 1*t;
  F.input("B")(2,0)   += 2*t;
  F.evaluate();
  cout << "solution (fd) = " << (F.output("X")-x_unpreturbed)/t << endl;

  // Jacobian
  Function J = F.jacobian("B","X");  
  J.init();
  J.setInput(val,"A");
  J.setInput(rhs,"B");
  J.evaluate();
  cout << "solution (dx/db) = " << J.output() << endl;
  DMatrix J_analytic = inv(J.input("A"));
  if(tr) J_analytic = J_analytic.T();
  cout << "analytic solution (dx/db) = " << J_analytic << endl;
}
