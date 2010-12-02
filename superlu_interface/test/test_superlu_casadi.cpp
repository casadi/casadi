/** 
File superlu.c from the SuperLU example collection
Joel Andersson, K.U. Leuven, 2010
*/

#include "superlu_interface/superlu.hpp"
#include "external_packages/superlu_4_1/SRC/slu_ddefs.h"

using namespace CasADi;

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
  SuperLU linear_solver(nrow,ncol,rowind,col);
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
  
  linear_solver.evaluate();
}
