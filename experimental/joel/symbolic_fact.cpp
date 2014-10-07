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


#include "core/std_vector_tools.hpp"
#include "core/function/sx_function.hpp"
#include "core/function/mx_function.hpp"
#include "core/function/external_function.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/mx/mx_tools.hpp"
#include "core/matrix/crs_sparsity_internal.hpp"
#include "core/matrix/sparsity_tools.hpp"
#include "interfaces/csparse/csparse.hpp"
extern "C"{
#include "external_packages/CSparse/Include/cs.h"
}

#include <ctime>

using namespace std;
using namespace casadi;

int main(){
  
  // Size
  int nrow = 50;
  int ncol = 50;
  
  // Density
  double dens = 0.12;
  
  // Number of nonzero elements
  int nnz = std::max(std::max(nrow,ncol),int(nrow*ncol*dens));
  
  // Seed
  double seed = 0.331;
  
  // Generate "random" (but deterministic) sparsity pattern
  vector<int> row(nnz), col(nnz);
  
  // Add diagonal
  for(int k=0; k<nrow; ++k){
    row[k] = k;
    col[k] = k;
  }
  
  for(int k=std::max(nrow,ncol); k<nnz; ++k){
    // Generate a random number between 0 and 1
    seed = 100*fabs(sin(seed*123 + 3));
    seed = seed - floor(seed);

    // Generate a random number between 0 and 300
    row[k] = std::min(int(seed*nrow),nrow-1);

    // Generate a random number between 0 and 1
    seed = 100*fabs(sin(seed*123 + 3));
    seed = seed - int(seed);

    // Generate a random number between 0 and 300
    col[k] = std::min(int(seed*ncol),ncol-1);
  }
  
  // Create sparsity pattern
  CRSSparsity C1 = sp_triplet(nrow, ncol, row, col);

  // This is the bcsstk01 matrix from CSparse
  int nrow2 = 48;
  int ncol2 = 48;
  int row2[224] = {0,4,5,6,10,18,24,29,1,3,5,7,9,19,23,25,2,3,4,8,20,22,26,27,3,7,9,21,26,27,4,6,10,20,22,28,5,11,19,23,24,29,6,10,11,12,30,35,7,9,11,13,17,31,8,9,10,14,16,32,33,9,15,32,33,10,14,16,34,11,13,17,30,
               35,12,16,17,18,22,36,41,42,46,47,13,14,15,17,19,21,37,43,44,45,14,15,16,20,38,39,43,44,45,15,19,21,38,39,43,44,45,16,17,18,22,40,42,46,47,17,23,36,41,42,46,47,18,22,23,42,47,19,21,23,43,20,21,22,
               44,45,21,44,45,22,46,23,42,47,24,28,29,30,34,25,27,31,33,26,27,32,27,31,33,28,30,34,29,35,30,34,35,36,31,33,35,37,41,32,33,34,38,40,33,39,34,38,40,35,37,41,36,40,41,42,46,37,39,41,43,45,38,39,40,
               44,39,43,45,40,42,46,41,47,42,46,47,43,44,45,44,45,45,46,47,47};
  int col2[224] = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,
               13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,18,18,18,18,18,19,19,19,19,20,20,20,20,20,21,21,21,22,22,23,23,23,24,24,24,
               24,24,25,25,25,25,26,26,26,27,27,27,28,28,28,29,29,30,30,30,30,31,31,31,31,31,32,32,32,32,32,33,33,34,34,34,35,35,35,36,36,36,36,36,37,37,37,37,37,38,38,38,38,39,39,39,40,40,40,41,41,42,42,42,43,
               43,43,44,44,45,46,46,47};
  double val2[224] = {2.83226851852e+06,1.0e+06,2.08333333333e+06,-3.33333333333e+03,1.0e+06,-2.8e+06,-2.89351851852e+04,2.08333333333e+06,1.63544753086e+06,-2.0e+06,5.55555555555e+06,-6.66666666667e+03,-2.0e+06,
                 -3.08641975309e+04,5.55555555555e+06,-1.59791666667e+06,1.72436728395e+06,-2.08333333333e+06,-2.77777777778e+06,-1.68e+06,-1.54320987654e+04,-2.77777777778e+06,-2.89351851852e+04,-2.08333333333e+06,1.00333333333e+09,2.0e+06,4.0e+08,-3.33333333333e+06,
                  2.08333333333e+06,1.0e+08,1.06750000000e+09,-1.0e+06,2.0e+08,2.77777777778e+06,3.33333333333e+08,-8.33333333333e+05,1.53533333333e+09,-2.0e+06,-5.55555555555e+06,6.66666666667e+08,
                  -2.08333333333e+06,1.0e+08,2.83226851852e+06,-1.0e+06,2.08333333333e+06,-2.8e+06,-2.89351851852e+04,2.08333333333e+06,1.63544753086e+06,2.0e+06,5.55555555555e+06,-3.08641975309e+04,5.55555555555e+06,-1.59791666667e+06,1.72436728395e+06,-2.08333333333e+06,-2.77777777778e+06,
                  -1.54320987654e+04,-2.77777777778e+06,-2.89351851852e+04,-2.08333333333e+06,1.00333333333e+09,-3.33333333333e+06,2.08333333333e+06,1.0e+08,1.06750000000e+09,2.77777777778e+06,3.33333333333e+08,
                  -8.33333333333e+05,1.53533333333e+09,-5.55555555555e+06,6.66666666667e+08,-2.08333333333e+06,1.0e+08,2.83609946950e+06,-2.14928529451e+06,2.35916180402e+06,-3.33333333333e+03,-1.0e+06,-2.89351851852e+04,2.08333333333e+06,-3.83095098171e+03,-1.14928529451e+06,
                  2.75828470683e+05,1.76741074446e+06,5.17922131816e+05,4.29857058902e+06,-5.55555555555e+06,-6.66666666667e+03,2.0e+06,-1.59791666667e+06,-1.31963213599e+05,-5.17922131816e+05,2.29857058902e+06,
                  3.89003806848e+06,-2.63499027470e+06,2.77777777778e+06,-1.68e+06,-2.89351851852e+04,-2.08333333333e+06,-5.17922131816e+05,-2.16567078453e+06,-5.51656941367e+05,1.97572063531e+09,-2.0e+06,4.0e+08,2.08333333333e+06,1.0e+08,-2.29857058902e+06,5.51656941366e+05,
                  4.86193650990e+08,1.52734651547e+09,-1.09779731332e+08,1.0e+06,2.0e+08,-8.33333333333e+05,1.14928529451e+06,2.29724661236e+08,-5.57173510779e+07,1.56411143711e+09,-2.0e+06,-2.08333333333e+06,1.0e+08,-2.75828470683e+05,-5.57173510779e+07,1.09411960038e+07,
                  2.83226851852e+06,1.0e+06,2.08333333333e+06,-2.89351851852e+04,2.08333333333e+06,1.63544753086e+06,-2.0e+06,-5.55555555555e+06,-1.59791666667e+06,1.72436728395e+06,-2.08333333333e+06,2.77777777778e+06,
                  -2.89351851852e+04,-2.08333333333e+06,1.00333333333e+09,2.08333333333e+06,1.0e+08,1.06750000000e+09,-8.33333333333e+05,1.53533333333e+09,-2.08333333333e+06,1.0e+08,6.08796296296e+04,1.25e+06,
                  4.16666666667e+05,-4.16666666667e+03,1.25e+06,3.37291666667e+06,-2.5e+06,-8.33333333333e+03,-2.5e+06,2.41171296296e+06,-4.16666666667e+05,-2.35500000000e+06,1.5e+09,2.5e+06,5.0e+08,5.01833333333e+08,-1.25e+06,2.5e+08,5.02500000000e+08,-2.5e+06,3.98587962963e+06,
                  -1.25e+06,4.16666666667e+05,-3.92500000000e+06,3.41149691358e+06,2.5e+06,6.94444444444e+06,-3.85802469136e+04,6.94444444445e+06,2.43100308642e+06,-4.16666666667e+05,-3.47222222222e+06,-1.92901234568e+04,-3.47222222222e+06,1.50416666667e+09,-4.16666666667e+06,
                  1.33516666667e+09,3.47222222222e+06,4.16666666667e+08,2.16916666667e+09,-6.94444444444e+06,8.33333333333e+08,3.98587962963e+06,-1.25e+06,4.16666666667e+05,-4.16666666667e+03,-1.25e+06,3.41149691358e+06,
                  2.5e+06,-6.94444444445e+06,-8.33333333333e+03,2.5e+06,2.43100308642e+06,-4.16666666667e+05,3.47222222222e+06,-2.35500000000e+06,1.50416666667e+09,-2.5e+06,5.0e+08,1.33516666667e+09,1.25e+06,2.5e+08,
                  2.16916666667e+09,-2.5e+06,6.47105806113e+04,2.39928529451e+06,1.40838195984e+05,3.50487988027e+06,5.17922131816e+05,-4.79857058902e+06,4.57738374749e+06,1.34990274700e+05,2.47238730198e+09,9.61679848804e+08,-1.09779731332e+08,5.31278103775e+08};

  vector<int> mapping;
  CRSSparsity C2 = sp_triplet(nrow2, ncol2, vector<int>(row2,row2+224), vector<int>(col2,col2+224),mapping);
  vector<double> valv(mapping.size());
  for(int k=0; k<valv.size(); ++k){
    valv[k] = val2[mapping[k]];
  }

  // Create a matrix
  DMatrix D2(C2,valv);
//  D2.printSparse();
    
  // Right hand side
  vector<double> b2(48,1.0);
    
  // Create a CSparse instance
  CSparse L2(C2);
  L2.init();
  L2.setInput(D2,0);
  L2.setInput(b2,1);
  L2.evaluate();
  
  
  // Now test a nonsymmetric matrix
  int rowind3[11] = {0,1,2,4,4,5,6,6,6,8,10};
  int col3[10] = {0,1,0,1,2,4,2,3,3,4};
  int nrow3 = 10;
  int ncol3 = 5;
  CRSSparsity S3(nrow3,ncol3,vector<int>(col3,col3+10),vector<int>(rowind3,rowind3+11));
  
  IMatrix(S3,1).printDense();

  
  cs AT_;
  AT_.nzmax = S3.size();  // maximum number of entries 
  AT_.m = S3.size2(); // number of rows
  AT_.n = S3.size1(); // number of columns
  AT_.p = &S3.rowindRef().front(); // column pointers (size n+1) or col indices (size nzmax)
  AT_.i = &S3.colRef().front(); // row indices, size nzmax
  AT_.x = 0; // row indices, size nzmax
  AT_.nz = -1; // of entries in triplet matrix, -1 for compressed-col 


  int order = 3;
  int qr = 1;
  css *S_ = cs_sqr(order, &AT_, qr);

  std::vector<int> pinv;
  std::vector<int> q;
  std::vector<int> parent; 
  std::vector<int> cp;
  std::vector<int> leftmost;
  int m2;
  double lnz;
  double unz;
  S3->prefactorize(order, qr, pinv, q, parent, cp, leftmost, m2, lnz, unz);

  cout << "pinv" << endl;
  cout << pinv << endl;
  if(S_->pinv!=0)
    cout << vector<int>(S_->pinv, S_->pinv + pinv.size()) << endl;
  cout << endl;
    
  cout << "q" << endl;
  cout << q << endl;
  if(S_->q!=0)
    cout << vector<int>(S_->q, S_->q + q.size()) << endl;
  cout << endl;

  cout << "parent" << endl;
  cout << parent << endl;
  if(S_->parent!=0)
    cout << vector<int>(S_->parent, S_->parent + parent.size()) << endl;
  cout << endl;

  cout << "cp" << endl;
  cout << cp << endl;
  if(S_->cp!=0)
    cout << vector<int>(S_->cp, S_->cp + cp.size()) << endl;
  cout << endl;
  
  cout << "leftmost" << endl;
  cout << leftmost << endl;
  if(S_->leftmost!=0)
    cout << vector<int>(S_->leftmost, S_->leftmost + leftmost.size()) << endl;
  cout << endl;
  
  cs_sfree(S_);
  
  int dmseed = 0;
  csd *perm = cs_dmperm (&AT_, dmseed);
    
  // Save to BLT structure // NOTE: swapping row<>col due to row/column major
  vector<int> rowperm3(perm->q, perm->q + S3.size1());
  vector<int> colperm3(perm->p, perm->p + S3.size2());
  int nb3 = perm->nb;
  vector<int> rowblock3(perm->s, perm->s + nb3 + 1);
  vector<int> colblock3(perm->r, perm->r + nb3 + 1);
  vector<int> coarse_rowblock3(perm->cc, perm->cc+5);
  vector<int> coarse_colblock3(perm->rr, perm->rr+5);
  
  cout << "csparse" << endl;
  cout << rowperm3 << endl;
  cout << colperm3 << endl;
  cout << rowblock3 << endl;
  cout << colblock3 << endl;
  cout << coarse_rowblock3 << endl;
  cout << coarse_colblock3 << endl;
  
  
  S3.dulmageMendelsohn(rowperm3, colperm3, rowblock3, colblock3, coarse_rowblock3, coarse_colblock3, dmseed);

  cout << "casadi" << endl;
  cout << rowperm3 << endl;
  cout << colperm3 << endl;
  cout << rowblock3 << endl;
  cout << colblock3 << endl;
  cout << coarse_rowblock3 << endl;
  cout << coarse_colblock3 << endl;

  
  // Free allocated memory and return
  cs_dfree(perm);

  
  return 0;
}
