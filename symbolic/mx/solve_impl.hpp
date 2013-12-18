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

#ifndef SOLVE_IMPL_HPP
#define SOLVE_IMPL_HPP

#include "solve.hpp"
#include "../matrix/matrix_tools.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include "../fx/fx_internal.hpp"

using namespace std;

namespace CasADi{

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A, const LinearSolver& linear_solver) : linear_solver_(linear_solver){
    casadi_assert_message(r.size2() == A.size1(),"Solve::Solve: dimension mismatch.");
    setDependencies(r,A);
    setSparsity(r.sparsity());
  }

  template<bool Tr>
  void Solve<Tr>::printPart(std::ostream &stream, int part) const{
    if(part==0){
      stream << "(";
    } else if(part==1){
      stream << "/";
    } else {
      if(Tr) stream << "'";
      stream << ")";
    }
  }

  template<bool Tr>
  void Solve<Tr>::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    
    // Factorize the matrix
    linear_solver_.setInput(*input[1],LINSOL_A);
    linear_solver_.prepare();
    
    // Solve for nondifferentiated output
    if(input[0]!=output[0]){
      copy(input[0]->begin(),input[0]->end(),output[0]->begin());
    }
    linear_solver_.solve(getPtr(output[0]->data()),output[0]->size1(),Tr);

    // Forward sensitivities
    for(int d=0; d<nfwd; ++d){
      if(fwdSeed[d][0]!=fwdSens[d][0]){
        copy(fwdSeed[d][0]->begin(),fwdSeed[d][0]->end(),fwdSens[d][0]->begin());
      }
      transform(fwdSens[d][0]->begin(),fwdSens[d][0]->end(),fwdSens[d][0]->begin(),std::negate<double>());
      if(Tr){
        DMatrix::mul_no_alloc_nt(*output[0],*fwdSeed[d][1],*fwdSens[d][0]);
      } else {
        DMatrix::mul_no_alloc_nn(*output[0],*fwdSeed[d][1],*fwdSens[d][0]);
      }
      transform(fwdSens[d][0]->begin(),fwdSens[d][0]->end(),fwdSens[d][0]->begin(),std::negate<double>());
      linear_solver_.solve(getPtr(fwdSens[d][0]->data()),output[0]->size1(),Tr);      
    }

    // Adjoint sensitivities
    for(int d=0; d<nadj; ++d){

      // Solve transposed
      transform(adjSeed[d][0]->begin(),adjSeed[d][0]->end(),adjSeed[d][0]->begin(),std::negate<double>());
      linear_solver_.solve(getPtr(adjSeed[d][0]->data()),output[0]->size1(),!Tr);

      // Propagate to A
      if(!Tr){
        DMatrix::mul_no_alloc_tn(*output[0],*adjSeed[d][0],*adjSens[d][1]);
      } else {
        DMatrix::mul_no_alloc_tn(*adjSeed[d][0],*output[0],*adjSens[d][1]);
      }

      // Propagate to B
      if(adjSeed[d][0]==adjSens[d][0]){
        transform(adjSens[d][0]->begin(),adjSens[d][0]->end(),adjSens[d][0]->begin(),std::negate<double>());
      } else {
        transform(adjSens[d][0]->begin(),adjSens[d][0]->end(),adjSeed[d][0]->begin(),adjSens[d][0]->begin(),std::minus<double>());
        fill(adjSeed[d][0]->begin(),adjSeed[d][0]->end(),0);
      }
    }
  }

  template<bool Tr>
  void Solve<Tr>::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
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
        X = A->getSolve(B,Tr,linear_solver_);
      }
    }

    // Forward sensitivities, collect the right hand sides
    vector<int> rhs_ind;
    vector<MX> rhs;
    vector<int> row_offset(1,0);
    for(int d=0; d<nfwd; ++d){
      const MX& B_hat = *fwdSeed[d][0];
      const MX& A_hat = *fwdSeed[d][1];
      
      // Get right hand side
      MX rhs_d;
      if(Tr){
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
      rhs = vertsplit(A->getSolve(vertcat(rhs),Tr,linear_solver_),row_offset);
    
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
      rhs = vertsplit(A->getSolve(vertcat(rhs),!Tr,linear_solver_),row_offset);
    
      for(int i=0; i<rhs.size(); ++i){
        int d = rhs_ind[i];

        // Propagate to A
        if(!Tr){
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
  
  template<bool Tr>
  void Solve<Tr>::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd){
#if 1

    // Sparsity of the rhs
    const CRSSparsity& sp = input[0]->sparsity();
    const vector<int>& rowind = sp.rowind();
    int nrhs = sp.size1();
    int nnz = input[1]->size();

    // Get pointers to data
    bvec_t* B_ptr = reinterpret_cast<bvec_t*>(input[0]->ptr());
    bvec_t* A_ptr = reinterpret_cast<bvec_t*>(input[1]->ptr());
    bvec_t* X_ptr = reinterpret_cast<bvec_t*>(output[0]->ptr());

    if(fwd){
    
      // Get dependencies of all A elements
      bvec_t A_dep = 0;
      for(int i=0; i<nnz; ++i){
        A_dep |= *A_ptr++;
      }

      // One right hand side at a time
      for(int i=0; i<nrhs; ++i){

        // Add dependencies on B
        bvec_t AB_dep = A_dep;
        for(int k=rowind[i]; k<rowind[i+1]; ++k){
          AB_dep |= *B_ptr++;
        }

        // Propagate to X
        for(int k=rowind[i]; k<rowind[i+1]; ++k){
          *X_ptr++ = AB_dep;
        }
      }
    } else {

      // Dependencies of all X
      bvec_t X_dep_all = 0;

      // One right hand side at a time
      for(int i=0; i<nrhs; ++i){

        // Everything that depends on X
        bvec_t X_dep = 0;
        for(int k=rowind[i]; k<rowind[i+1]; ++k){
          X_dep |= *X_ptr;
          *X_ptr++ = 0;
        }

        // Propagate to B
        for(int k=rowind[i]; k<rowind[i+1]; ++k){
          *B_ptr++ |= X_dep;
        }

        // Collect dependencies on all X
        X_dep_all |= X_dep;
      }

      // Propagate to A
      for(int i=0; i<nnz; ++i){
        *A_ptr++ |= X_dep_all;
      }
    }


#else

    // Sparsities
    const CRSSparsity& r_sp = input[0]->sparsity();
    const CRSSparsity& A_sp = input[1]->sparsity();
    const vector<int>& A_rowind = A_sp.rowind();
    const vector<int>& A_col = A_sp.col();
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

        // Copy B_ptr to a temporary vector and clear X_ptr (NOTE: normally, B_ptr == X_ptr)
        copy(B_ptr,B_ptr+n,tmp_ptr);
        std::fill(X_ptr,X_ptr+n,0);

        // Add A_hat contribution to tmp
        for(int i=0; i<n; ++i){ // Loop over the rows of A
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
            int j = A_col[k]; // get the column
            if(Tr){
              tmp_ptr[i] |= A_ptr[k];
            } else {
              tmp_ptr[j] |= A_ptr[k];
            }
          }
        }
                
        // Propagate dependencies to X_hat
        for(int i=0; i<n; ++i){ // Loop over the rows of A
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
            int j = A_col[k]; // get the column            
            if(Tr){
              X_ptr[i] |= tmp_ptr[j];
            } else {
              X_ptr[j] |= tmp_ptr[i];
            }
          }
        }
        
        // Propagate interdependencies in X_hat
        std::fill(tmp_ptr,tmp_ptr+n,0);
        for(int iter=0; iter<n; ++iter){ // n represents the worst case
          bool no_change = true; // is there any change at all in this iteration
          for(int i=0; i<n; ++i){ // Loop over the rows of A
            for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
              int j = A_col[k]; // get the column
              if(Tr){
                no_change &= X_ptr[i] == tmp_ptr[i];
                X_ptr[i] |= tmp_ptr[i];
                tmp_ptr[i] |= X_ptr[i];
              } else {
                no_change &= X_ptr[j] == tmp_ptr[j];
                X_ptr[j] |= tmp_ptr[j];
                tmp_ptr[j] |= X_ptr[j];
              }
            }
          }

          // Stop iterating if reached a steady state
          if(no_change) break;
        }
      
      } else { // adjoint

        // Propagate dependencies from X_ptr
        std::fill(tmp_ptr,tmp_ptr+n,0);
        for(int i=0; i<n; ++i){ // Loop over the rows of A
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
            int j = A_col[k]; // get the column            
            if(Tr){
              tmp_ptr[j] |= X_ptr[i];
            } else {
              tmp_ptr[i] |= X_ptr[j];
            }
          }
        }

        // Propagate interdependencies, use X_ptr as temporary
        std::fill(X_ptr,X_ptr+n,0);
        for(int iter=0; iter<n; ++iter){ // n represents the worst case
          bool no_change = true; // is there any change at all in this iteration
          for(int i=0; i<n; ++i){ // Loop over the rows of A
            for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
              int j = A_col[k]; // get the column
              if(Tr){
                no_change &= tmp_ptr[i] == X_ptr[i];
                tmp_ptr[i] |= X_ptr[i];
                X_ptr[i] |= tmp_ptr[i];
              } else {
                no_change &= tmp_ptr[j] == X_ptr[j];
                tmp_ptr[j] |= X_ptr[j];
                X_ptr[j] |= tmp_ptr[j];
              }
            }
          }

          // Stop iterating if reached a steady state
          if(no_change) break;
        }

        // Clear seeds (again)
        std::fill(X_ptr,X_ptr+n,0);

        // Propagate to A_ptr
        for(int i=0; i<n; ++i){ // Loop over the rows of A
          for(int k=A_rowind[i]; k<A_rowind[i+1]; ++k){ // loop over the nonzeros
            int j = A_col[k]; // get the column
            if(Tr){
              A_ptr[k] |= tmp_ptr[i];
            } else {
              A_ptr[k] |= tmp_ptr[j];
            }
          }
        }

        // Propagate to B_ptr
        for(int i=0; i<n; ++i){ // Loop over the rows of A
          B_ptr[i] |= tmp_ptr[i];
        }
      }

      // Continue to the next right-hand-side
      B_ptr += n;
      X_ptr += n;      
    }

#endif
  }

  template<bool Tr>
  void Solve<Tr>::deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied) {
    MXNode::deepCopyMembers(already_copied);
    linear_solver_ = deepcopy(linear_solver_, already_copied);
  }

} // namespace CasADi

#endif // SOLVE_IMPL_HPP

