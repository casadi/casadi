/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "superlu_internal.hpp"

using namespace std;
namespace CasADi{
  
SuperLUInternal::SuperLUInternal(int nrow, int ncol, const vector<int>& rowind, const vector<int>& col, int nrhs) : 
  nrow_(nrow), ncol_(ncol), rowind_(rowind), col_(col), nrhs_(nrhs){
    
  // Allocate space for inputs
  input_.resize(2);
  input_[0].setSize(nrow,ncol);
  input_[0].setSparsityCRS(rowind, col);
  
  input_[1].setSize(nrow,nrhs); // right hand side
  
  // Allocate space for outputs
  output_.resize(1);
  output_[0].setSize(ncol,nrhs);
  
  // not initialized
  is_init = false;
}

SuperLUInternal::~SuperLUInternal(){
  if(is_init){
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    SUPERLU_FREE (rhs);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    StatFree(&stat);
  }
}

void SuperLUInternal::init(){
  // Call the init method of the base class
  FXNode::init();
  
  // Make sure not already initialized
  if(is_init) throw CasadiException("SuperLU: Already initialized");
  
  // number of non-zero elements
  int nnz = col_.size();

  // Allocate SuperLU data structures
  a = doubleMalloc(nnz);
  asub = intMalloc(nnz);
  xa = intMalloc(ncol_+1);
  rhs = doubleMalloc(nrow_ * nrhs_);
  perm_r = intMalloc(nrow_);
  perm_c = intMalloc(ncol_);
  
  if(a==0 || asub==0 || xa==0 || rhs==0 || perm_r==0 || perm_c==0) 
    throw CasadiException("SuperLU: Malloc failed.");

  // Copy the values
  copy(col_.begin(),col_.end(),asub);
  copy(rowind_.begin(),rowind_.end(),xa);

  // Set to initialized
  is_init = true;
  
}

void SuperLUInternal::evaluate(int fsens_order, int asens_order){
  // number of non-zero elements
  int nnz = col_.size();

  // Copy the non-zero entries
  const vector<double>& val = input(0).data();
  copy(val.begin(),val.end(),a);

  // Copy the right hand side
  const vector<double>& b = input(1).data();
  copy(b.begin(),b.end(),rhs);

  // Create matrix A in the format expected by SuperLU.
  dCreate_CompCol_Matrix(&A, nrow_, ncol_, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&B, nrow_, nrhs_, rhs, nrow_, SLU_DN, SLU_D, SLU_GE);

  // Set the default input options
  set_default_options(&options);
  options.ColPerm = NATURAL;

  // Initialize the statistics variables
  StatInit(&stat);

  // Solve the linear system.
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
  dPrint_CompCol_Matrix(const_cast<char*>("A"), &A);
  dPrint_CompCol_Matrix(const_cast<char*>("U"), &U);
  dPrint_SuperNode_Matrix(const_cast<char*>("L"), &L);
  print_int_vec(const_cast<char*>("\nperm_r"), nrow_, perm_r);

}
  
  
} // namespace CasADi

  


