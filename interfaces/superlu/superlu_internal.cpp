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

#include "superlu_internal.hpp"

using namespace std;
namespace CasADi{
  
SuperLUInternal::SuperLUInternal(const CRSSparsity& sparsity, int nrhs)  : LinearSolverInternal(sparsity), nrhs_(nrhs){
  
  // Add options
  addOption("equil", OT_BOOLEAN, true); // Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).
  addOption("colperm", OT_STRING, "colamd"); // Specifies how to permute the columns of the matrix for sparsity preservation.
  addOption("iterrefine", OT_STRING, "norefine"); 
  addOption("diagpivotthresh",OT_REAL,1.0);
  addOption("symmetricmode",OT_BOOLEAN,false);
  addOption("pivotgrowth",OT_BOOLEAN,false);
  addOption("conditionnumber",OT_BOOLEAN,false);
  addOption("rowperm",OT_STRING,"largediag");
  addOption("printstat",OT_BOOLEAN,true);
  addOption("user_work", OT_BOOLEAN, false); // keep work in memory
      
  // not initialized
  is_init_ = false;
  
  // Work vector
  work_ = 0;
  lwork_ = 0;
}

SuperLUInternal::SuperLUInternal(const SuperLUInternal& linsol) : LinearSolverInternal(linsol){
  // not initialized
  is_init_ = false;
  
  // Work vector
  work_ = 0;
  lwork_ = 0;
}

SuperLUInternal::~SuperLUInternal(){
  unInit();
  if(work_) free(work_);
}

void SuperLUInternal::unInit(){
  if(is_init_){
//    Destroy_CompCol_Matrix(&A_); // Not allowed since we allocate all data in stl arrays and pass as pointers
//    Destroy_SuperMatrix_Store(&B_); // Not allowed since we allocate all data in stl arrays and pass as pointers
    Destroy_SuperNode_Matrix(&L_);
    Destroy_CompCol_Matrix(&U_);
    StatFree(&stat_);
    is_init_=false;
  }
}

void SuperLUInternal::init(){
  // Call the init method of the base class
  LinearSolverInternal::init();

  // Reuse work
  user_work_ = getOption("user_work").toInt();

  // Free allocated memory, if any
  unInit();
  
  // Make sure not already initialized
  if(is_init_) throw CasadiException("SuperLU: Already initialized");

  // Make sure that the matrix is is square (relax later)
  if(nrow()!=ncol()) throw CasadiException("SuperLU: Not square");
  
  // Allocate SuperLU data structures
  a_.resize(nnz());
  rhs_.resize(nrow() * nrhs_);
  perm_r_.resize(nrow());
  perm_c_.resize(ncol());

  
  // Create matrices A and B in the format expected by SuperLU.
  dCreate_CompCol_Matrix(&A_, ncol(), nrow(), nnz(), &a_[0], const_cast<int*>(&col()[0]), const_cast<int*>(&rowind()[0]), SLU_NC, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&B_, nrow(), nrhs_, &rhs_[0], nrow(), SLU_DN, SLU_D, SLU_GE);

  // Initialize the statistics variables
  StatInit(&stat_);

  // Set the default input options
  set_default_options(&options_);

  // Read column permutation
  GenericType colperm = getOption("colperm");
  if(colperm=="natural")            options_.ColPerm = ::NATURAL;
  else if(colperm=="mmd_ata")       options_.ColPerm = ::MMD_ATA;
  else if(colperm=="mmd_at_plus_a") options_.ColPerm = ::MMD_AT_PLUS_A;
  else if(colperm=="colamd")        options_.ColPerm = ::COLAMD;
  else if(colperm=="my_permc")      options_.ColPerm = ::MY_PERMC;
  else throw CasadiException("SuperLU: Unknown column permutation: " + colperm.toString());

  // Read transpose
  options_.Trans = transpose_ ? ::NOTRANS : ::TRANS; // swap due to row-major/col-major

  // Iterataive refinement
  GenericType iterrefine = getOption("iterrefine");
  if(iterrefine=="norefine" || iterrefine=="no")  options_.IterRefine = ::NOREFINE; // user guide is inconsistent, allow both possibilties
  else if(iterrefine=="single")                   options_.IterRefine = ::SINGLE;
  else if(iterrefine=="double")                   options_.IterRefine = ::DOUBLE;
  else if(iterrefine=="extra")                    options_.IterRefine = ::EXTRA;
  else throw CasadiException("SuperLU: Unknown iterative refinement: " + iterrefine.toString());

  // Specifies the threshold used for a diagonal entry to be an acceptable pivot.
  options_.DiagPivotThresh = getOption("diagpivotthresh").toDouble();
  
  // Specifies whether to use the symmetric mode. Symmetric mode gives preference to diagonal pivots, and uses an (AT + A)-based column permutation algorithm.
  options_.SymmetricMode = getOption("symmetricmode").toInt() ? YES : NO;

  // Specifies whether to compute the reciprocal pivot growth.
  options_.PivotGrowth = getOption("pivotgrowth").toInt() ? YES : NO;

  // Specifies whether to compute the reciprocal condition number.
  options_.ConditionNumber = getOption("conditionnumber").toInt() ? YES : NO;
  
  // Specifies whether to permute the rows of the original matrix.
  GenericType rowperm = getOption("rowperm");
  if(rowperm=="no" || rowperm=="norowperm")    options_.RowPerm = ::NOROWPERM; 
  else if(rowperm=="largediag")                options_.RowPerm = ::LargeDiag;
  else if(rowperm=="my_permr")                 options_.RowPerm = ::MY_PERMR;
  else throw CasadiException("SuperLU: Unknown row permutation: " + rowperm.toString());

  // Specifies which version of modified ILU to use.
  options_.PrintStat = getOption("printstat").toInt() ? YES : NO;
  
  // Elimination tree
  etree_.resize(A_.ncol);
  
  // Set to initialized
  is_init_ = true;
  
  // Has the routine been called once
  called_once_ = false;

  // Internal work
  work_ = 0;

  // length of work vector
  lwork_ = 0; 
}

void SuperLUInternal::prepare(){
  prepared_ = false;

  // Copy the non-zero entries
  const vector<double>& val = input(0);
  copy(val.begin(),val.end(),a_.begin());

  SuperMatrix AC;
  get_perm_c(options_.ColPerm, &A_, &perm_c_[0]);
  sp_preorder(&options_, &A_, &perm_c_[0], &etree_[0], &AC);

  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2); // no of columns in a relaxed snodes

  // If no work array, estimate the needed size
if (user_work_){
  if(!lwork_){
    int sz = 0;
    lwork_ = -1;
    dgstrf(&options_, &AC, relax, panel_size, &etree_[0], work_, lwork_, &perm_c_[0], &perm_r_[0], &L_, &U_, &stat_, &sz);
    lwork_ = sz; // - A->ncol;
    if(work_) free(work_);
    work_ = malloc(lwork_);
  }
} // user_work_    
    
  // Compute the LU factorization of A
  do{
    info_ = 0;
    dgstrf(&options_, &AC, relax, panel_size, &etree_[0], work_, lwork_, &perm_c_[0], &perm_r_[0], &L_, &U_, &stat_, &info_);

    if(info_<0){
      stringstream ss;
      ss << "SuperLU: The " << (-info_) << "-th argument had an illegal value" << endl;
      throw CasadiException(ss.str());
    } else if(info_>0){
      if(info_<=A_.ncol){
        stringstream ss;
        ss << "SuperLU: U(" << info_ << "," << info_ << ") is exactly zero. "
        "The factorization has been completed, but the factor U is exactly singular, "
        "and division by zero will occur if it is used to solve a system of equations.";
        throw CasadiException(ss.str());
      } else {
        if (user_work_){
          // Allocate more memory and repeat
          cout << "Allocating more memory" << endl;
          lwork_ *= 2;
          work_ = realloc(work_,lwork_);
        } else { // user_work_
          stringstream ss;
          ss << "SuperLU: Allocation failed after " << (info_-A_.ncol) << " bytes allocated";
          throw CasadiException(ss.str());
        } // user_work_    
      }
    }
  } while(info_!=0);
  
  // Destroy temporary memory
  Destroy_CompCol_Permuted(&AC);

  called_once_ = true;
  prepared_ = true;
}
  
void SuperLUInternal::solve(double* x, int nrhs){
  casadi_assert(nrhs_==nrhs);
  
  // Copy the right hand side
  copy(x,x+rhs_.size(),rhs_.begin());
    
  // Solve the system A*X=B, overwriting B with X. 
  info_ = 0;
  dgstrs(options_.Trans, &L_, &U_, &perm_c_[0], &perm_r_[0], &B_, &stat_, &info_);
  if(info_ != 0) throw CasadiException("dgstrs failed");

  // Copy the result
  copy(rhs_.begin(),rhs_.end(),x);
}

SuperLUInternal* SuperLUInternal::clone() const{
  return new SuperLUInternal(*this);
}

  
} // namespace CasADi
