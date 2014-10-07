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


#include "linear_solver_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

INPUTSCHEME(LinsolInput)
OUTPUTSCHEME(LinsolOutput)

using namespace std;
namespace casadi {

  LinearSolverInternal::LinearSolverInternal(const Sparsity& sparsity, int nrhs) {
    // Make sure arguments are consistent
    casadi_assert(!sparsity.isNull());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "LinearSolverInternal::init: the matrix must be square but got "
                          << sparsity.dimString());
    casadi_assert_message(!sparsity.isSingular(),
                          "LinearSolverInternal::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << rank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    std::vector<int> coarse_rowblock, coarse_colblock;
    sparsity.dulmageMendelsohn(rowperm_, colperm_, rowblock_, colblock_,
                               coarse_rowblock, coarse_colblock);

    // Allocate inputs
    setNumInputs(LINSOL_NUM_IN);
    input(LINSOL_A) = DMatrix(sparsity);
    input(LINSOL_B) = DMatrix::zeros(sparsity.size2(), nrhs);

    // Allocate outputs
    setNumOutputs(LINSOL_NUM_OUT);
    output(LINSOL_X) = input(LINSOL_B);

    input_.scheme = SCHEME_LinsolInput;
    output_.scheme = SCHEME_LinsolOutput;
  }

  void LinearSolverInternal::init() {
    // Call the base class initializer
    FunctionInternal::init();

    // Not prepared
    prepared_ = false;
  }

  LinearSolverInternal::~LinearSolverInternal() {
  }

  void LinearSolverInternal::evaluate() {
    /*  Factorization fact;
        if (called_once) {
        // Check if any element has changed
        bool any_change = false;
        const vector<double>& val = input(0).data();
        for (int i=0; i<val.size(); ++i) {
        if (val[i] != a[i]) {
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
    if (!prepared_)
      throw CasadiException("LinearSolverInternal::evaluate: Preparation failed");

    // Solve the factorized system
    solve(false);
  }

  void LinearSolverInternal::solve(bool transpose) {
    // Get input and output vector
    const vector<double>& b = input(LINSOL_B).data();
    vector<double>& x = output(LINSOL_X).data();
    int nrhs = input(LINSOL_B).size2();

    // Copy input to output
    copy(b.begin(), b.end(), x.begin());

    // Solve the factorized system in-place
    solve(getPtr(x), nrhs, transpose);
  }

  void LinearSolverInternal::evaluateMXGen(const MXPtrV& input, MXPtrV& output,
                                           const MXPtrVV& fwdSeed, MXPtrVV& fwdSens,
                                           const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                                           bool output_given, bool tr) {
    int nfwd = fwdSens.size();
    int nadj = adjSeed.size();
    const MX& B = *input[0];
    const MX& A = *input[1];
    MX& X = *output[0];

    // Nondifferentiated output
    if (!output_given) {
      if (B.isZero()) {
        X = MX::sparse(B.shape());
      } else {
        X = solve(A, B, tr);
      }
    }

    // Forward sensitivities, collect the right hand sides
    std::vector<int> rhs_ind;
    std::vector<MX> rhs;
    std::vector<int> col_offset(1, 0);
    for (int d=0; d<nfwd; ++d) {
      const MX& B_hat = *fwdSeed[d][0];
      const MX& A_hat = *fwdSeed[d][1];

      // Get right hand side
      MX rhs_d;
      if (tr) {
        rhs_d = B_hat - mul(A_hat.T(), X);
      } else {
        rhs_d = B_hat - mul(A_hat, X);
      }

      // Simplifiy if zero
      if (rhs_d.isZero()) {
        *fwdSens[d][0] = MX::sparse(rhs_d.shape());
      } else {
        rhs.push_back(rhs_d);
        rhs_ind.push_back(d);
        col_offset.push_back(col_offset.back()+rhs_d.size2());
      }
    }

    if (!rhs.empty()) {
      // Solve for all directions at once
      rhs = horzsplit(solve(A, horzcat(rhs), tr), col_offset);

      // Save result
      for (int i=0; i<rhs.size(); ++i) {
        *fwdSens[rhs_ind[i]][0] = rhs[i];
      }
    }

    // Adjoint sensitivities, collect right hand sides
    rhs.resize(0);
    rhs_ind.resize(0);
    col_offset.resize(1);
    for (int d=0; d<nadj; ++d) {
      MX& X_bar = *adjSeed[d][0];

      // Simplifiy if zero
      if (X_bar.isZero()) {
        if (adjSeed[d][0]!=adjSens[d][0]) {
          *adjSens[d][0] = X_bar;
          X_bar = MX();
        }
      } else {
        rhs.push_back(X_bar);
        rhs_ind.push_back(d);
        col_offset.push_back(col_offset.back()+X_bar.size2());

        // Delete seed
        X_bar = MX();
      }
    }

    if (!rhs.empty()) {
      // Solve for all directions at once
      rhs = horzsplit(solve(A, horzcat(rhs), !tr), col_offset);

      for (int i=0; i<rhs.size(); ++i) {
        int d = rhs_ind[i];

        // Propagate to A
        if (!tr) {
          adjSens[d][1]->addToSum(-mul(rhs[i], X.T(), A.sparsity()));
        } else {
          adjSens[d][1]->addToSum(-mul(X, rhs[i].T(), A.sparsity()));
        }

        // Propagate to B
        if (adjSeed[d][0]==adjSens[d][0]) {
          *adjSens[d][0] = rhs[i];
        } else {
          adjSens[d][0]->addToSum(rhs[i]);
        }
      }
    }
  }

  void LinearSolverInternal::propagateSparsityGen(DMatrixPtrV& input, DMatrixPtrV& output,
                                                  std::vector<int>& itmp, std::vector<double>& rtmp,
                                                  bool fwd, bool transpose) {

    // Sparsities
    const Sparsity& r_sp = input[0]->sparsity();
    const Sparsity& A_sp = input[1]->sparsity();
    const std::vector<int>& A_colind = A_sp.colind();
    const std::vector<int>& A_row = A_sp.row();
    int nrhs = r_sp.size2();
    int n = r_sp.size1();
    //    int nnz = A_sp.size();

    // Get pointers to data
    bvec_t* B_ptr = reinterpret_cast<bvec_t*>(input[0]->ptr());
    bvec_t* A_ptr = reinterpret_cast<bvec_t*>(input[1]->ptr());
    bvec_t* X_ptr = reinterpret_cast<bvec_t*>(output[0]->ptr());
    bvec_t* tmp_ptr = reinterpret_cast<bvec_t*>(getPtr(rtmp));

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {

      if (fwd) {

        // Copy B_ptr to a temporary vector
        copy(B_ptr, B_ptr+n, tmp_ptr);

        // Add A_hat contribution to tmp
        for (int cc=0; cc<n; ++cc) {
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr = A_row[k];
            tmp_ptr[transpose ? cc : rr] |= A_ptr[k];
          }
        }

        // Propagate to X_ptr
        std::fill(X_ptr, X_ptr+n, 0);
        spSolve(X_ptr, tmp_ptr, transpose);

      } else { // adjoint

        // Solve transposed
        std::fill(tmp_ptr, tmp_ptr+n, 0);
        spSolve(tmp_ptr, B_ptr, !transpose);

        // Clear seeds
        std::fill(B_ptr, B_ptr+n, 0);

        // Propagate to X_ptr
        for (int i=0; i<n; ++i) {
          X_ptr[i] |= tmp_ptr[i];
        }

        // Propagate to A_ptr
        for (int cc=0; cc<n; ++cc) {
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr = A_row[k];
            A_ptr[k] |= tmp_ptr[transpose ? cc : rr];
          }
        }
      }

      // Continue to the next right-hand-side
      B_ptr += n;
      X_ptr += n;
    }
  }

  void LinearSolverInternal::spSolve(DMatrix& X, const DMatrix& B, bool transpose) const {
    bvec_t* X_bvec = reinterpret_cast<bvec_t*>(X.ptr());
    const bvec_t* B_bvec = reinterpret_cast<const bvec_t*>(B.ptr());
    spSolve(X_bvec, B_bvec, transpose);
  }

  void LinearSolverInternal::spSolve(bvec_t* X, const bvec_t* B, bool transpose) const {

    const Sparsity& A_sp = input(LINSOL_A).sparsity();
    const std::vector<int>& A_colind = A_sp.colind();
    const std::vector<int>& A_row = A_sp.row();
    int nb = rowblock_.size()-1; // number of blocks

    if (!transpose) {
      for (int b=0; b<nb; ++b) { // loop over the blocks forward

        // Get dependencies from all right-hand-sides in the block ...
        bvec_t block_dep = 0;
        for (int el=rowblock_[b]; el<rowblock_[b+1]; ++el) {
          int rr = rowperm_[el];
          block_dep |= B[rr];
        }

        // ... as well as all other variables in the block
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];
          block_dep |= X[cc];
        }

        // Propagate ...
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];

          // ... to all variables in the block ...
          X[cc] |= block_dep;

          // ... as well as to other variables which depends on variables in the block
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr=A_row[k];
            X[rr] |= block_dep;
          }
        }
      }

    } else { // transpose
      for (int b=nb-1; b>=0; --b) { // loop over the blocks backward

        // Get dependencies ...
        bvec_t block_dep = 0;
        for (int el=colblock_[b]; el<colblock_[b+1]; ++el) {
          int cc = colperm_[el];

          // .. from all right-hand-sides in the block ...
          block_dep |= B[cc];

          // ... as well as from all depending variables ...
          for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
            int rr=A_row[k];
            block_dep |= X[rr];
          }
        }

        // Propagate to all variables in the block
        for (int el=rowblock_[b]; el<rowblock_[b+1]; ++el) {
          int rr = rowperm_[el];
          X[rr] |= block_dep;
        }
      }
    }
  }

  void LinearSolverInternal::evaluateDGen(const DMatrixPtrV& input, DMatrixPtrV& output, bool tr) {

    // Factorize the matrix
    setInput(*input[1], LINSOL_A);
    prepare();

    // Solve for nondifferentiated output
    if (input[0]!=output[0]) {
      copy(input[0]->begin(), input[0]->end(), output[0]->begin());
    }
    solve(getPtr(output[0]->data()), output[0]->size2(), tr);
  }

  void LinearSolverInternal::evaluateSXGen(const SXPtrV& input, SXPtrV& output, bool tr) {
    casadi_error("LinearSolverInternal::evaluateSXGen not defined for class "
                 << typeid(*this).name());
  }

  MX LinearSolverInternal::solve(const MX& A, const MX& B, bool transpose) {
    return A->getSolve(B, transpose, shared_from_this<LinearSolver>());
  }

  void LinearSolverInternal::solve(double* x, int nrhs, bool transpose) {
    casadi_error("LinearSolverInternal::solve not defined for class " << typeid(*this).name());
  }

  std::map<std::string, LinearSolverInternal::Plugin> LinearSolverInternal::solvers_;

  const std::string LinearSolverInternal::infix_ = "linearsolver";

  void LinearSolverInternal::solveL(double* x, int nrhs, bool transpose) {
    casadi_error("LinearSolverInternal::solveL not defined for class "
                 << typeid(*this).name());
  }

  Sparsity LinearSolverInternal::getFactorizationSparsity(bool transpose) const {
    casadi_error("LinearSolverInternal::getFactorizationSparsity not defined for class "
                 << typeid(*this).name());
    return Sparsity();
  }

  DMatrix LinearSolverInternal::getFactorization(bool transpose) const {
    casadi_error("LinearSolverInternal::getFactorization not defined for class "
                 << typeid(*this).name());
    return DMatrix();
  }

} // namespace casadi

