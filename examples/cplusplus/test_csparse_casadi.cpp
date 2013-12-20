/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "symbolic/casadi.hpp"
#include "interfaces/csparse/csparse.hpp"

using namespace CasADi;
using namespace std;

int main(int argc, char *argv[])
{

  /* Initialize matrix A. */
  int nrow = 5, ncol = 5;
  int nnz = 12;
  
  vector<int> rowind(nrow+1);
  vector<int> col(nnz);
  
  // Sparsity pattern
  col[0] = 0; col[1] = 1; col[2] = 4; col[3] = 1;
  col[4] = 2; col[5] = 4; col[6] = 0; col[7] = 2;
  col[8] = 0; col[9] = 3; col[10]= 3; col[11]= 4;
  rowind[0] = 0; rowind[1] = 3; rowind[2] = 6; rowind[3] = 8; rowind[4] = 10; rowind[5] = 12;
  CRSSparsity spA(nrow,ncol,col,rowind);
  
  // Create a solver instance
  CSparse linear_solver(spA);
    
  // Initialize
  linear_solver.init();

  // Pass Non-zero elements
  double   s, u, p, e, r, l;
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;

  vector<double> val(nnz);
  val[0] = s; val[1] = l; val[2] = l; val[3] = u; val[4] = l; val[5] = l;
  val[6] = u; val[7] = p; val[8] = u; val[9] = e; val[10]= u; val[11]= r;
  
  // Right hand side
  vector<double> rhs(nrow,1.0);
  
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
  MX A = msym("A",spA);
  MX B = msym("B",1,nrow);
  MX X = linear_solver.solve(A,B,tr);
  MXFunction F(linsolIn("A",A,"B",B),linsolOut("X",X));
  F.init();

  // Solve nondifferentiated
  F.setInput(val,"A");
  F.setInput(rhs,"B");
  F.fwdSeed("A")(2,3)   = 1;
  F.fwdSeed("B")(0,2)   = 2;
  F.adjSeed("X")(0,3)   = 1;
  F.evaluate(1,1);
  
  // Print the solution
  cout << "solution (embedded) = " << F.output("X") << endl;
  cout << "forward sensitivities = " << F.fwdSens("X") << endl;
  cout << "adjoint sensitivities (A) = " << F.adjSens("A") << endl;
  cout << "adjoint sensitivities (B) = " << F.adjSens("B") << endl;
  
  // Preturb the linear solver
  double t = 0.01;
  DMatrix x_unpreturbed = F.output("X");
  F.input("A")(2,3)   += 1*t;
  F.input("B")(0,2)   += 2*t;
  F.evaluate();
  cout << "solution (fd) = " << (F.output("X")-x_unpreturbed)/t << endl;

  // Jacobian
  FX J = F.jacobian("B","X");  
  J.init();
  J.setInput(val,"A");
  J.setInput(rhs,"B");
  J.evaluate();
  cout << "solution (dx/db) = " << J.output() << endl;
  DMatrix J_analytic = inv(J.input("A"));
  if(!tr) J_analytic = trans(J_analytic);
  cout << "analytic solution (dx/db) = " << J_analytic << endl;
}
