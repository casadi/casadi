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
  
  // Add options
  addOption("equil", OT_BOOLEAN, true); // Specifies whether to equilibrate the system (scale A’s rows and columns to have unit norm).
  addOption("colperm", OT_STRING, "colamd"); // Specifies how to permute the columns of the matrix for sparsity preservation.
  addOption("trans", OT_BOOLEAN, false); // transpose A
  addOption("iterrefine", OT_STRING, "norefine"); 
  addOption("diagpivotthresh",OT_REAL,1.0);
  addOption("symmetricmode",OT_BOOLEAN,false);
  addOption("pivotgrowth",OT_BOOLEAN,false);
  addOption("conditionnumber",OT_BOOLEAN,false);
  addOption("rowperm",OT_STRING,"largediag");
  addOption("printstat",OT_BOOLEAN,true);
    
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
    Destroy_CompRow_Matrix(&A);
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

  // Create matrices A and B in the format expected by SuperLU.
  dCreate_CompRow_Matrix(&A, nrow_, ncol_, nnz, a, asub, xa, SLU_NR, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&B, nrow_, nrhs_, rhs, nrow_, SLU_DN, SLU_D, SLU_GE);
 
  // Initialize the statistics variables
  StatInit(&stat);

  // Set the default input options
  set_default_options(&options);
  
  // Read column permutation
  Option colperm = getOption("colperm");
  if(colperm=="natural")            options.ColPerm = ::NATURAL;
  else if(colperm=="mmd_ata")       options.ColPerm = ::MMD_ATA;
  else if(colperm=="mmd_at_plus_a") options.ColPerm = ::MMD_AT_PLUS_A;
  else if(colperm=="colamd")        options.ColPerm = ::COLAMD;
  else if(colperm=="my_permc")      options.ColPerm = ::MY_PERMC;
  else throw CasadiException("SuperLU: Unknown column permutation: " + colperm.toString());
  
  // Read transpose
  options.Trans = getOption("trans").toInt() ? ::TRANS : ::NOTRANS;
  
  // Iterataive refinement
  Option iterrefine = getOption("iterrefine");
  if(iterrefine=="norefine" || iterrefine=="no")  options.IterRefine = ::NOREFINE; // user guide is inconsistent, allow both possibilties
  else if(iterrefine=="single")                   options.IterRefine = ::SINGLE;
  else if(iterrefine=="double")                   options.IterRefine = ::DOUBLE;
  else if(iterrefine=="extra")                    options.IterRefine = ::EXTRA;
  else throw CasadiException("SuperLU: Unknown iterative refinement: " + iterrefine.toString());

  // Specifies the threshold used for a diagonal entry to be an acceptable pivot.
  options.DiagPivotThresh = getOption("diagpivotthresh").toDouble();
  
  // Specifies whether to use the symmetric mode. Symmetric mode gives preference to diagonal pivots, and uses an (AT + A)-based column permutation algorithm.
  options.SymmetricMode = getOption("symmetricmode").toInt() ? YES : NO;

  // Specifies whether to compute the reciprocal pivot growth.
  options.PivotGrowth = getOption("pivotgrowth").toInt() ? YES : NO;

  // Specifies whether to compute the reciprocal condition number.
  options.ConditionNumber = getOption("conditionnumber").toInt() ? YES : NO;
  
  // Specifies whether to permute the rows of the original matrix.
  Option rowperm = getOption("rowperm");
  if(rowperm=="no" || rowperm=="norowperm")    options.RowPerm = ::NOROWPERM; 
  else if(rowperm=="largediag")                options.RowPerm = ::LargeDiag;
  else if(rowperm=="my_permr")                 options.RowPerm = ::MY_PERMR;
  else throw CasadiException("SuperLU: Unknown row permutation: " + rowperm.toString());

  // Specifies which version of modified ILU to use.
  options.PrintStat = getOption("printstat").toInt() ? YES : NO;
  
  // Set to initialized
  is_init = true;
  
}

void SuperLUInternal::evaluate(int fsens_order, int asens_order){
  solve(DOFACT);
}

void SuperLUInternal::solve(Factorization fact){
  // Copy the non-zero entries
  const vector<double>& val = input(0).data();
  copy(val.begin(),val.end(),a);

  // Copy the right hand side
  const vector<double>& b = input(1).data();
  copy(b.begin(),b.end(),rhs);

  // Choose factorization
  if(fact==DOFACT)                        options.Fact = ::DOFACT;
  else if(fact==SAMEPATTERN)              options.Fact = ::SamePattern;
  else if(fact==SAMEPATTERN_SAMEROWPERM)  options.Fact = ::SamePattern_SameRowPerm;
  else if(fact==FACTORED)                 options.Fact = ::FACTORED;

  // Solve the linear system
  dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

  // Copy the result
  vector<double>& res = output(0).data();
  copy(rhs,rhs+res.size(),res.begin());
  
  dPrint_CompCol_Matrix(const_cast<char*>("A"), &A);
  dPrint_CompCol_Matrix(const_cast<char*>("U"), &U);
  dPrint_SuperNode_Matrix(const_cast<char*>("L"), &L);
  print_int_vec(const_cast<char*>("\nperm_r"), nrow_, perm_r);
  
}

  
} // namespace CasADi

  


