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


#include "linsol_internal.hpp"
#include "../std_vector_tools.hpp"
#include "../mx/mx_node.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  int linsol_n_in() {
    return LINSOL_NUM_IN;
  }

  int linsol_n_out() {
    return LINSOL_NUM_OUT;
  }

  vector<string> linsol_in() {
    vector<string> ret(linsol_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=linsol_in(i);
    return ret;
  }

  vector<string> linsol_out() {
    vector<string> ret(linsol_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=linsol_out(i);
    return ret;
  }

  string linsol_in(int ind) {
    switch (static_cast<LinsolInput>(ind)) {
    case LINSOL_A:     return "A";
    case LINSOL_B:      return "B";
    case LINSOL_NUM_IN: break;
    }
    return string();
  }

  string linsol_out(int ind) {
    switch (static_cast<LinsolOutput>(ind)) {
    case LINSOL_X:     return "X";
    case LINSOL_NUM_OUT: break;
    }
    return string();
  }

  bool has_linsol(const string& name) {
    return LinsolInternal::has_plugin(name);
  }

  void load_linsol(const string& name) {
    LinsolInternal::load_plugin(name);
  }

  string doc_linsol(const string& name) {
    return LinsolInternal::getPlugin(name).doc;
  }

  Function linsol_new(const std::string& name, const std::string& solver,
                  const Sparsity& sp, int nrhs, const Dict& opts) {
    Linsol F(name + "_linsol", solver, sp, opts);
    MX A = MX::sym("A", sp);
    MX b = MX::sym("b", sp.size2(), nrhs);
    MX x = F.solve(A, b);
    return Function(name, {A, b}, {x}, {"A", "B"}, {"X"});
  }

  Function linsol(const std::string& name, const std::string& solver,
                  const Sparsity& sp, int nrhs, const Dict& opts) {
    casadi_assert(nrhs==0);
    Function ret;
    if (solver=="none") {
      ret.assignNode(new LinsolInternal(name, sp, nrhs));
    } else {
      ret.assignNode(LinsolInternal::getPlugin(solver).creator(name, sp, nrhs));
    }
    ret->construct(opts);
    return ret;
  }

  DM Function::linsol_solve(const DM& A, const DM& B, bool tr) {
    // Factorize
    linsol_factorize(A.ptr());

    // Solve
    DM x = densify(B);
    linsol_solve(x.ptr(), x.size2());
    return x;
  }

  MX Function::linsol_solve(const MX& A, const MX& B, bool tr) {
    return (*this)->linsol_solve(A, B, tr);
  }

  void Function::linsol_solveL(double* x, int nrhs, bool tr, int mem) const {
    (*this)->linsol_solveL(memory(mem), x, nrhs, tr);
  }

  void Function::linsol_factorize(const double* A, int mem) const {
    (*this)->linsol_factorize(memory(mem), A);
  }

  void Function::linsol_solve(double* x, int nrhs, bool tr, int mem) const {
    (*this)->linsol_solve(memory(mem), x, nrhs, tr);
  }

  Sparsity Function::linsol_cholesky_sparsity(bool tr, int mem) const {
    return (*this)->linsol_cholesky_sparsity(memory(mem), tr);
  }

  DM Function::linsol_cholesky(bool tr, int mem) const {
    return (*this)->linsol_cholesky(memory(mem), tr);
  }

  LinsolInternal::LinsolInternal(const std::string& name, const Sparsity& sparsity, int nrhs)
    : FunctionInternal(name), sparsity_(sparsity), nrhs_(nrhs) {

    // Make sure arguments are consistent
    casadi_assert(!sparsity.is_null());
    casadi_assert_message(sparsity.size2()==sparsity.size1(),
                          "LinsolInternal::init: the matrix must be square but got "
                          << sparsity.dim());
    casadi_assert_message(!sparsity.is_singular(),
                          "LinsolInternal::init: singularity - the matrix is structurally "
                          "rank-deficient. sprank(J)=" << sprank(sparsity)
                          << " (in stead of "<< sparsity.size2() << ")");

    // Calculate the Dulmage-Mendelsohn decomposition
    btf_ = sparsity.btf();

    // Number of equations
    neq_ = sparsity.size2();
  }

  LinsolInternal::~LinsolInternal() {
  }

  Sparsity LinsolInternal::get_sparsity_in(int i) {
    switch (static_cast<LinsolInput>(i)) {
    case LINSOL_A:
      return sparsity_;
    case LINSOL_B:
      return Sparsity::dense(neq_, nrhs_);
    case LINSOL_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity LinsolInternal::get_sparsity_out(int i) {
    switch (static_cast<LinsolOutput>(i)) {
    case LINSOL_X:
      return Sparsity::dense(neq_, nrhs_);
    case LINSOL_NUM_OUT: break;
    }
    return Sparsity();
  }

  void LinsolInternal::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

  }

  void LinsolInternal::eval(void* mem, const double** arg, double** res,
                    int* iw, double* w) const {
    // Get inputs and outputs
    const double *A = arg[LINSOL_A];
    const double *b = arg[LINSOL_B];
    arg += LINSOL_NUM_IN;
    double *x = res[LINSOL_X];
    res += LINSOL_NUM_OUT;

    // If output not requested, nothing to do
    if (!x) return;

    // A zero linear system would be singular
    if (A==0) {
      casadi_fill(x, neq_*nrhs_, numeric_limits<double>::quiet_NaN());
      return;
    }

    // If right hand side is zero, solution is trivially zero (if well-defined)
    if (!b) {
      casadi_fill(x, neq_*nrhs_, 0.);
      return;
    }

    // Setup memory object
    setup(mem, arg, res, iw, w);

    // Factorize the linear system
    linsol_factorize(mem, A);

    // Solve the factorized system
    casadi_copy(b, neq_*nrhs_, x);
    linsol_solve(mem, x, nrhs_, false);
  }

  void LinsolInternal::
  linsol_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& fsens, bool tr) {
    // Number of derivatives
    int nfwd = fseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nfwd);
    std::vector<int> col_offset(nfwd+1, 0);
    for (int d=0; d<nfwd; ++d) {
      const MX& B_hat = fseed[d][0];
      const MX& A_hat = fseed[d][1];
      rhs[d] = tr ? B_hat - mtimes(A_hat.T(), X) : B_hat - mtimes(A_hat, X);
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_solve(A, horzcat(rhs), tr), col_offset);

    // Fetch result
    fsens.resize(nfwd);
    for (int d=0; d<nfwd; ++d) {
      fsens[d].resize(1);
      fsens[d][0] = rhs[d];
    }
  }

  void LinsolInternal::
  linsol_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                 const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& asens, bool tr) {
    // Number of derivatives
    int nadj = aseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nadj);
    std::vector<int> col_offset(nadj+1, 0);
    for (int d=0; d<nadj; ++d) {
      rhs[d] = aseed[d][0];
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(linsol_solve(A, horzcat(rhs), !tr), col_offset);

    // Collect sensitivities
    asens.resize(nadj);
    for (int d=0; d<nadj; ++d) {
      asens[d].resize(2);

      // Propagate to A
      MX a;
      if (!tr) {
        a = -mac(rhs[d], X.T(), MX::zeros(A.sparsity()));
      } else {
        a = -mac(X, rhs[d].T(), MX::zeros(A.sparsity()));
      }
      if (asens[d][1].is_empty(true)) {
        asens[d][1] = a;
      } else {
        asens[d][1] += a;
      }

      // Propagate to B
      if (asens[d][0].is_empty(true)) {
        asens[d][0] = rhs[d];
      } else {
        asens[d][0] += rhs[d];
      }
    }
  }

  void LinsolInternal::
  linsol_sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_in(LINSOL_A);
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    const bvec_t *B=arg[0], *A = arg[1];
    bvec_t* X = res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {
      // Copy B to a temporary vector
      copy(B, B+n, tmp);

      // Add A_hat contribution to tmp
      for (int cc=0; cc<n; ++cc) {
        for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          int rr = A_row[k];
          tmp[tr ? cc : rr] |= A[k];
        }
      }

      // Propagate to X
      std::fill(X, X+n, 0);
      A_sp.spsolve(btf_, X, tmp, tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void LinsolInternal::
  linsol_sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
               bool tr, int nrhs) {
    // Sparsities
    const Sparsity& A_sp = sparsity_in(LINSOL_A);
    const int* A_colind = A_sp.colind();
    const int* A_row = A_sp.row();
    int n = A_sp.size1();

    // Get pointers to data
    bvec_t *B=arg[0], *A=arg[1], *X=res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (int r=0; r<nrhs; ++r) {
      // Solve transposed
      std::fill(tmp, tmp+n, 0);
      A_sp.spsolve(btf_, tmp, X, !tr);

      // Clear seeds
      std::fill(X, X+n, 0);

      // Propagate to B
      for (int i=0; i<n; ++i) B[i] |= tmp[i];

      // Propagate to A
      for (int cc=0; cc<n; ++cc) {
        for (int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          int rr = A_row[k];
          A[k] |= tmp[tr ? cc : rr];
        }
      }

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
  }

  void LinsolInternal::linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                             bool tr, int nrhs) {
    casadi_error("LinsolInternal::eval_sxLinsol not defined for class "
                 << typeid(*this).name());
  }

  MX LinsolInternal::linsol_solve(const MX& A, const MX& B, bool tr) {
    return A->getSolve(B, tr, shared_from_this<Function>());
  }

  std::map<std::string, LinsolInternal::Plugin> LinsolInternal::solvers_;

  const std::string LinsolInternal::infix_ = "linsol";

} // namespace casadi
