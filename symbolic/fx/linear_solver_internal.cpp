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

#include "linear_solver_internal.hpp"
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../mx/mx_node.hpp"

INPUTSCHEME(LinsolInput)
OUTPUTSCHEME(LinsolOutput)

using namespace std;
namespace CasADi{

  LinearSolverInternal::LinearSolverInternal(const CRSSparsity& sparsity, int nrhs){
    // Make sure arguments are consistent
    casadi_assert(!sparsity.isNull());
    casadi_assert_message(sparsity.size1()==sparsity.size2(),"LinearSolverInternal::init: the matrix must be square but got " << sparsity.dimString());  
    casadi_assert_message(!isSingular(sparsity),"LinearSolverInternal::init: singularity - the matrix is structurally rank-deficient. sprank(J)=" << rank(sparsity) << " (in stead of "<< sparsity.size1() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    std::vector<int> coarse_rowblock, coarse_colblock;
    sparsity.dulmageMendelsohn(rowperm_, colperm_, rowblock_, colblock_, coarse_rowblock, coarse_colblock);

    // Allocate inputs
    setNumInputs(LINSOL_NUM_IN);
    input(LINSOL_A) = DMatrix(sparsity);
    input(LINSOL_B) = DMatrix(nrhs,sparsity.size1(),0);
  
    // Allocate outputs
    setNumOutputs(LINSOL_NUM_OUT);
    output(LINSOL_X) = input(LINSOL_B);

    input_.scheme = SCHEME_LinsolInput;
    output_.scheme = SCHEME_LinsolOutput;
  }

  void LinearSolverInternal::init(){
    // Call the base class initializer
    FXInternal::init();
    
    // Not prepared
    prepared_ = false;
  }

  LinearSolverInternal::~LinearSolverInternal(){
  }
 
  void LinearSolverInternal::evaluate(){
    /*  Factorization fact;
        if(called_once){
        // Check if any element has changed
        bool any_change = false;
        const vector<double>& val = input(0).data();
        for(int i=0; i<val.size(); ++i){
        if(val[i] != a[i]){
        any_change = true;
        break;
        }
        }
    
        // Reuse factored matrix if matrix hasn't changed
        fact = any_change ? SAMEPATTERN : FACTORED;
        } else {
        fact = DOFACT;
        called_once = true;
        }*/
  
    // Call the solve routine
    prepare();
  
    // Make sure preparation successful
    if(!prepared_) 
      throw CasadiException("LinearSolverInternal::evaluate: Preparation failed");
  
    // Solve the factorized system
    solve(false);
  }
 
  void LinearSolverInternal::solve(bool transpose){
    // Get input and output vector
    const vector<double>& b = input(LINSOL_B).data();
    vector<double>& x = output(LINSOL_X).data();
    int nrhs = input(LINSOL_B).size1();

    // Copy input to output
    copy(b.begin(),b.end(),x.begin());
  
    // Solve the factorized system in-place
    solve(getPtr(x),nrhs,transpose);
  }

  void LinearSolverInternal::evaluateMXGen(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given, bool tr){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    const MX& B = *input[0];
    const MX& A = *input[1];
    MX& X = *output[0];

    // Nondifferentiated output
    if(!output_given){
      if(CasADi::isZero(B)){
        X = MX::sparse(B.shape());
      } else {
        X = solve(A,B,tr);
      }
    }

    // Forward sensitivities, collect the right hand sides
    std::vector<int> rhs_ind;
    std::vector<MX> rhs;
    std::vector<int> row_offset(1,0);
    for(int d=0; d<nfwd; ++d){
      const MX& B_hat = *fwdSeed[d][0];
      const MX& A_hat = *fwdSeed[d][1];
      
      // Get right hand side
      MX rhs_d;
      if(tr){
        rhs_d = B_hat - mul(X,trans(A_hat));
      } else {
        rhs_d = B_hat - mul(X,A_hat);
      }
      
      // Simplifiy if zero
      if(CasADi::isZero(rhs_d)){
        *fwdSens[d][0] = MX::sparse(rhs_d.shape());
      } else {
        rhs.push_back(rhs_d);
        rhs_ind.push_back(d);
        row_offset.push_back(row_offset.back()+rhs_d.size1());
      }
    }
    
    if(!rhs.empty()){
      // Solve for all directions at once
      rhs = vertsplit(solve(A,vertcat(rhs),tr),row_offset);
    
      // Save result
      for(int i=0; i<rhs.size(); ++i){
        *fwdSens[rhs_ind[i]][0] = rhs[i];
      }
    }

    // Adjoint sensitivities, collect right hand sides
    rhs.resize(0);
    rhs_ind.resize(0);
    row_offset.resize(1);
    for(int d=0; d<nadj; ++d){
      MX& X_bar = *adjSeed[d][0];
      
      // Simplifiy if zero
      if(CasADi::isZero(X_bar)){
        if(adjSeed[d][0]!=adjSens[d][0]){
          *adjSens[d][0] = X_bar;
          X_bar = MX();
        }
      } else {
        rhs.push_back(X_bar);
        rhs_ind.push_back(d);
        row_offset.push_back(row_offset.back()+X_bar.size1());

        // Delete seed
        X_bar = MX();
      }
    }

    if(!rhs.empty()){
      // Solve for all directions at once
      rhs = vertsplit(solve(A,vertcat(rhs),!tr),row_offset);
    
      for(int i=0; i<rhs.size(); ++i){
        int d = rhs_ind[i];

        // Propagate to A
        if(!tr){
          *adjSens[d][1] -= mul(trans(X),rhs[i],A.sparsity());
        } else {
          *adjSens[d][1] -= mul(trans(rhs[i]),X,A.sparsity());
        }

        // Propagate to B
        if(adjSeed[d][0]==adjSens[d][0]){
          *adjSens[d][0] = rhs[i];
        } else {
          *adjSens[d][0] += rhs[i];
        }
      }
    }
  }

  void LinearSolverInternal::propagateSparsityGen(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd, bool transpose){

    // Sparsities
    const CRSSparsity& r_sp = input[0]->sparsity();
    const CRSSparsity& A_sp = input[1]->sparsity();
    const std::vector<int>& A_rowind = A_sp.rowind();
    const std::vector<int>& A_col = A_sp.col();
    int nrhs = r_sp.size1();
    int n = r_sp.size2();
    int nnz = A_sp.size();

    // Get pointers to data
    bvec_t* B_ptr = reinterpret_cast<bvec_t*>(input[0]->ptr());
    bvec_t* A_ptr = reinterpret_cast<bvec_t*>(input[1]->ptr());
    bvec_t* X_ptr = reinterpret_cast<bvec_t*>(output[0]->ptr());
    bvec_t* tmp_ptr = reinterpret_cast<bvec_t*>(getPtr(rtmp));

    // For all right-hand-sides
    for(int r=0; r<nrhs; ++r){

      if(fwd){

        // Copy B_ptr to a temporary vector
        copy(B_ptr,B_ptr+n,tmp_ptr);

        // Add A_hat contribution to tmp
        for(int i=0; i<n; ++i){
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){
            int j = A_col[k];
            tmp_ptr[transpose ? i : j] |= A_ptr[k];
          }
        }

        // Propagate to X_ptr
        std::fill(X_ptr,X_ptr+n,0);
        spSolve(X_ptr,tmp_ptr,transpose);
      
      } else { // adjoint
        
        // Solve transposed
        std::fill(tmp_ptr,tmp_ptr+n,0);
        spSolve(tmp_ptr,B_ptr,!transpose);
        
        // Clear seeds
        std::fill(B_ptr,B_ptr+n,0);

        // Propagate to X_ptr
        for(int i=0; i<n; ++i){
          X_ptr[i] |= tmp_ptr[i];
        }

        // Propagate to A_ptr
        for(int i=0; i<n; ++i){
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){
            int j = A_col[k];
            A_ptr[k] |= tmp_ptr[transpose ? i : j];
          }
        }
      }

      // Continue to the next right-hand-side
      B_ptr += n;
      X_ptr += n;      
    }
  }

  void LinearSolverInternal::spSolve(bvec_t* X, bvec_t* B, bool transpose) const{
    const CRSSparsity& A_sp = input(LINSOL_A).sparsity();
    const std::vector<int>& A_rowind = A_sp.rowind();
    const std::vector<int>& A_col = A_sp.col();

    if(transpose){
      int nb = colblock_.size()-1; // number of blocks
      for(int b=0; b<nb; ++b){ // loop over the blocks

        // Get dependencies ...
        bvec_t block_dep = 0;
        for(int el=rowblock_[b]; el<rowblock_[b+1]; ++el){
          // ... from all right-hand-sides in the block ...
          int i = rowperm_[el];
          block_dep |= B[i];

          // ... as well as from all dependent variables
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){
            int j=A_col[k];
            block_dep |= X[j];
          }
        }

        // Propagate to all variables in the block
        for(int el=colblock_[b]; el<colblock_[b+1]; ++el){
          int j = colperm_[el];
          X[j] = block_dep;
        }
      }
    } else {
      int nb = colblock_.size()-1; // number of blocks
      for(int b=nb-1; b>=0; --b){ // loop over the blocks
            
        // Get dependencies from all right-hand-sides
        bvec_t block_dep = 0;
        for(int el=colblock_[b]; el<colblock_[b+1]; ++el){
          int j = colperm_[el];
          block_dep |= B[j];
        }

        // Propagate ...
        for(int el=rowblock_[b]; el<rowblock_[b+1]; ++el){
          // ... to all variables in the block
          int i = rowperm_[el];
          X[i] = block_dep;

          // ... as well as to all other right-hand-sides
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){
            int j=A_col[k];
            B[j] |= block_dep;
          }
        }
      }
    }
  }

  void LinearSolverInternal::evaluateDGen(const DMatrixPtrV& input, DMatrixPtrV& output, bool tr){
    
    // Factorize the matrix
    setInput(*input[1],LINSOL_A);
    prepare();
    
    // Solve for nondifferentiated output
    if(input[0]!=output[0]){
      copy(input[0]->begin(),input[0]->end(),output[0]->begin());
    }
    solve(getPtr(output[0]->data()),output[0]->size1(),tr);
  }

  MX LinearSolverInternal::solve(const MX& A, const MX& B, bool transpose){
    return A->getSolve(B, transpose, shared_from_this<LinearSolver>());
  }
 
} // namespace CasADi
 
