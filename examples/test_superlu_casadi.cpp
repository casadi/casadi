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
File superlu.c from the SuperLU example collection
Joel Andersson, K.U. Leuven, 2010
*/

#include "casadi/stl_vector_tools.hpp"
#include "interfaces/superlu/superlu.hpp"

using namespace CasADi;
using namespace std;

main(int argc, char *argv[])
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
  
  // Create a solver instance
  SuperLU linear_solver(CRSSparsity(nrow,ncol,col,rowind));
  
  // Set options
  linear_solver.setOption("colperm", "natural");
  
  
  // Initialize
  linear_solver.init();

  // Pass Non-zero elements
  double   s, u, p, e, r, l;
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;

  vector<double> val(nnz);
  val[0] = s; val[1] = l; val[2] = l; val[3] = u; val[4] = l; val[5] = l;
  val[6] = u; val[7] = p; val[8] = u; val[9] = e; val[10]= u; val[11]= r;
  linear_solver.setInput(val,0);
  
  // Pass right hand side
  vector<double> rhs(nrow,1.0);
  linear_solver.setInput(rhs,1);
  
  // Solve
  linear_solver.evaluate();
  
  // Print the solution
  cout << "solution = " << linear_solver.output() << endl;
  
}
