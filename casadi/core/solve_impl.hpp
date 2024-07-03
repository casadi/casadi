/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_SOLVE_IMPL_HPP
#define CASADI_SOLVE_IMPL_HPP

#include "solve.hpp"
#include "linsol_internal.hpp"

namespace casadi {

  template<bool Tr>
  Solve<Tr>::Solve(const MX& r, const MX& A) {
    casadi_assert(r.size1() == A.size2(),
      "Solve::Solve: dimension mismatch. Got r " + r.dim() + " and A " + A.dim());
    set_dep(r, A);
    set_sparsity(r.sparsity());
  }

  template<bool Tr>
  std::string Solve<Tr>::disp(const std::vector<std::string>& arg) const {
    std::stringstream ss;
    ss << "(" << mod_prefix() << arg.at(1) << mod_suffix();
    if (Tr) ss << "'";
    ss << "\\" << arg.at(0) << ")";
    return ss.str();
  }

  template<bool Tr>
  LinsolCall<Tr>::LinsolCall(const MX& r, const MX& A, const Linsol& linear_solver) :
      Solve<Tr>(r, A), linsol_(linear_solver) {
  }

  template<bool Tr>
  int LinsolCall<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    scoped_checkout<Linsol> mem(linsol_);

    auto m = static_cast<LinsolMemory*>(linsol_->memory(mem));
    // Reset statistics
    for (auto&& s : m->fstats) s.second.reset();
    if (m->t_total) m->t_total->tic();

    if (linsol_.sfact(arg[1], mem)) return 1;
    if (linsol_.nfact(arg[1], mem)) return 1;
    if (linsol_.solve(arg[1], res[0], this->dep(0).size2(), Tr, mem)) return 1;

    linsol_->print_time(m->fstats);

    return 0;
  }

  template<bool Tr>
  int LinsolCall<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    linsol_->linsol_eval_sx(arg, res, iw, w, linsol_->memory(0), Tr, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  void Solve<Tr>::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    if (arg[0].is_zero()) {
      res[0] = MX(arg[0].size());
    } else {
      res[0] = solve(arg[1], arg[0], Tr);
    }
  }

  template<bool Tr>
  void Solve<Tr>::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    // Nondifferentiated inputs and outputs
    std::vector<MX> arg(this->n_dep());
    for (casadi_int i=0; i<arg.size(); ++i) arg[i] = this->dep(i);
    std::vector<MX> res(this->nout());
    for (casadi_int i=0; i<res.size(); ++i) res[i] = this->get_output(i);

    // Number of derivatives
    casadi_int nfwd = fseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nfwd);
    std::vector<casadi_int> col_offset(nfwd+1, 0);
    for (casadi_int d=0; d<nfwd; ++d) {
      const MX& B_hat = fseed[d][0];
      const MX& A_hat = fseed[d][1];
      rhs[d] = Tr ? B_hat - mtimes(A_hat.T(), X) : B_hat - mtimes(A_hat, X);
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(solve(A, horzcat(rhs), Tr), col_offset);

    // Fetch result
    fsens.resize(nfwd);
    for (casadi_int d=0; d<nfwd; ++d) {
      fsens[d].resize(1);
      fsens[d][0] = rhs[d];
    }
  }

  template<bool Tr>
  void Solve<Tr>::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    // Nondifferentiated inputs and outputs
    std::vector<MX> arg(this->n_dep());
    for (casadi_int i=0; i<arg.size(); ++i) arg[i] = this->dep(i);
    std::vector<MX> res(this->nout());
    for (casadi_int i=0; i<res.size(); ++i) res[i] = this->get_output(i);

    // Number of derivatives
    casadi_int nadj = aseed.size();
    const MX& A = arg[1];
    const MX& X = res[0];

    // Solve for all directions at once
    std::vector<MX> rhs(nadj);
    std::vector<casadi_int> col_offset(nadj+1, 0);
    for (casadi_int d=0; d<nadj; ++d) {
      rhs[d] = aseed[d][0];
      col_offset[d+1] = col_offset[d] + rhs[d].size2();
    }
    rhs = horzsplit(solve(A, horzcat(rhs), !Tr), col_offset);

    // Collect sensitivities
    asens.resize(nadj);
    for (casadi_int d=0; d<nadj; ++d) {
      asens[d].resize(2);

      // Propagate to A
      MX a;
      if (!Tr) {
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

  template<bool Tr>
  int Solve<Tr>::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Number of right-hand-sides
    casadi_int nrhs = dep(0).size2();

    // Sparsities
    const Sparsity& A_sp = this->A_sp();
    const casadi_int* A_colind = A_sp.colind();
    const casadi_int* A_row = A_sp.row();
    casadi_int n = A_sp.size1();

    // Get pointers to data
    const bvec_t *B=arg[0], *A = arg[1];
    bvec_t* X = res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (casadi_int r=0; r<nrhs; ++r) {
      // Copy B to a temporary vector
      std::copy(B, B+n, tmp);

      // Add A_hat contribution to tmp
      for (casadi_int cc=0; cc<n; ++cc) {
        for (casadi_int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          casadi_int rr = A_row[k];
          tmp[Tr ? cc : rr] |= A[k];
        }
      }

      // Propagate to X
      std::fill(X, X+n, 0);
      A_sp.spsolve(X, tmp, Tr);

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
    return 0;
  }

  template<bool Tr>
  int Solve<Tr>::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    // Number of right-hand-sides
    casadi_int nrhs = dep(0).size2();

    // Sparsities
    const Sparsity& A_sp = this->A_sp();
    const casadi_int* A_colind = A_sp.colind();
    const casadi_int* A_row = A_sp.row();
    casadi_int n = A_sp.size1();

    // Get pointers to data
    bvec_t *B=arg[0], *A=arg[1], *X=res[0];
    bvec_t* tmp = w;

    // For all right-hand-sides
    for (casadi_int r=0; r<nrhs; ++r) {
      // Solve transposed
      std::fill(tmp, tmp+n, 0);
      A_sp.spsolve(tmp, X, !Tr);

      // Clear seeds
      std::fill(X, X+n, 0);

      // Propagate to B
      for (casadi_int i=0; i<n; ++i) B[i] |= tmp[i];

      // Propagate to A
      for (casadi_int cc=0; cc<n; ++cc) {
        for (casadi_int k=A_colind[cc]; k<A_colind[cc+1]; ++k) {
          casadi_int rr = A_row[k];
          A[k] |= tmp[Tr ? cc : rr];
        }
      }

      // Continue to the next right-hand-side
      B += n;
      X += n;
    }
    return 0;
  }

  template<bool Tr>
  size_t LinsolCall<Tr>::sz_w() const {
    return this->sparsity().size1();
  }

  template<bool Tr>
  void LinsolCall<Tr>::generate(CodeGenerator& g,
                            const std::vector<casadi_int>& arg,
                            const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = this->dep(0).size2();

    // Array for x
    g.local("rr", "casadi_real", "*");
    g << "rr = " << g.work(res[0], this->nnz()) << ";\n";

    // Array for A
    g.local("ss", "casadi_real", "*");
    g << "ss = " << g.work(arg[1], this->dep(1).nnz()) << ";\n";

    // Copy b to x if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->nnz()), this->nnz(), "rr") << '\n';
    }
    // Solver specific codegen
    linsol_->generate(g, "ss", "rr", nrhs, Tr);
  }

  template<bool Tr>
  void Solve<Tr>::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
  }

  template<bool Tr>
  void Solve<Tr>::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Solve::Tr", Tr);
  }

  template<bool Tr>
  Solve<Tr>::Solve(DeserializingStream& s) : MXNode(s) {
  }

  template<bool Tr>
  MXNode* Solve<Tr>::deserialize(DeserializingStream& s) {
    bool tr;
    s.unpack("Solve::Tr", tr);
    casadi_error("Not implemented");
  }

  template<bool Tr>
  void LinsolCall<Tr>::serialize_body(SerializingStream& s) const {
    Solve<Tr>::serialize_body(s);
    s.pack("Solve::Linsol", linsol_);
  }

  template<bool Tr>
  void LinsolCall<Tr>::serialize_type(SerializingStream& s) const {
    Solve<Tr>::serialize_type(s);
  }

  template<bool Tr>
  LinsolCall<Tr>::LinsolCall(DeserializingStream& s) : Solve<Tr>(s) {
    s.unpack("Solve::Linsol", linsol_);
  }

  template<bool Tr>
  MXNode* LinsolCall<Tr>::deserialize(DeserializingStream& s) {
    bool tr;
    s.unpack("Solve::Tr", tr);

    if (tr) {
      return new LinsolCall<true>(s);
    } else {
      return new LinsolCall<false>(s);
    }
  }

  template<bool Tr>
  TriuSolve<Tr>::TriuSolve(const MX& r, const MX& A) : Solve<Tr>(r, A) {
  }

  template<bool Tr>
  int TriuSolve<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_triusolve(this->dep(1).sparsity(), arg[1], res[0], Tr, false, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  int TriuSolve<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_triusolve(this->dep(1).sparsity(), arg[1], res[0], Tr, false, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  TrilSolve<Tr>::TrilSolve(const MX& r, const MX& A) : Solve<Tr>(r, A) {
  }

  template<bool Tr>
  int TrilSolve<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_trilsolve(this->dep(1).sparsity(), arg[1], res[0], Tr, false, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  int TrilSolve<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_trilsolve(this->dep(1).sparsity(), arg[1], res[0], Tr, false, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  SolveUnity<Tr>::SolveUnity(const MX& r, const MX& A) : Solve<Tr>(r, A) {
  }

  template<bool Tr>
  const Sparsity& SolveUnity<Tr>::A_sp() const {
    // Create on first call
    if (A_sp_.is_null()) {
      const Sparsity& no_diag = this->dep(1).sparsity();
      A_sp_ = no_diag + Sparsity::diag(no_diag.size1());
    }
    // Return reference
    return A_sp_;
  }

  template<bool Tr>
  TriuSolveUnity<Tr>::TriuSolveUnity(const MX& r, const MX& A)
    : SolveUnity<Tr>(r, A) {
  }

  template<bool Tr>
  int TriuSolveUnity<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_triusolve(this->dep(1).sparsity(), arg[1], res[0], Tr, true, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  int TriuSolveUnity<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw,
      SXElem* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_triusolve(this->dep(1).sparsity(), arg[1], res[0], Tr, true, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  TrilSolveUnity<Tr>::TrilSolveUnity(const MX& r, const MX& A)
    : SolveUnity<Tr>(r, A) {
  }

  template<bool Tr>
  int TrilSolveUnity<Tr>::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_trilsolve(this->dep(1).sparsity(), arg[1], res[0], Tr, true, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  int TrilSolveUnity<Tr>::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw,
      SXElem* w) const {
    if (arg[0] != res[0]) std::copy(arg[0], arg[0] + this->dep(0).nnz(), res[0]);
    casadi_trilsolve(this->dep(1).sparsity(), arg[1], res[0], Tr, true, this->dep(0).size2());
    return 0;
  }

  template<bool Tr>
  void TriuSolve<Tr>::generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = this->dep(0).size2();
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->nnz()), this->nnz(), g.work(res[0], this->nnz())) << '\n';
    }
    // Perform sparse matrix multiplication
    g << g.triusolve(this->dep(1).sparsity(), g.work(arg[1], this->dep(1).nnz()),
      g.work(res[0], this->nnz()), Tr, false, nrhs) << '\n';
  }

  template<bool Tr>
  void TrilSolve<Tr>::generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = this->dep(0).size2();
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->nnz()), this->nnz(), g.work(res[0], this->nnz())) << '\n';
    }
    // Perform sparse matrix multiplication
    g << g.trilsolve(this->dep(1).sparsity(), g.work(arg[1], this->dep(1).nnz()),
      g.work(res[0], this->nnz()), Tr, false, nrhs) << '\n';
  }

  template<bool Tr>
  void TriuSolveUnity<Tr>::generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = this->dep(0).size2();
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->nnz()), this->nnz(), g.work(res[0], this->nnz())) << '\n';
    }
    // Perform sparse matrix multiplication
    g << g.triusolve(this->dep(1).sparsity(), g.work(arg[1], this->dep(1).nnz()),
      g.work(res[0], this->nnz()), Tr, true, nrhs) << '\n';
  }

  template<bool Tr>
  void TrilSolveUnity<Tr>::generate(CodeGenerator& g, const std::vector<casadi_int>& arg,
      const std::vector<casadi_int>& res) const {
    // Number of right-hand-sides
    casadi_int nrhs = this->dep(0).size2();
    // Copy first argument if not inplace
    if (arg[0]!=res[0]) {
      g << g.copy(g.work(arg[0], this->nnz()), this->nnz(), g.work(res[0], this->nnz())) << '\n';
    }
    // Perform sparse matrix multiplication
    g << g.trilsolve(this->dep(1).sparsity(), g.work(arg[1], this->dep(1).nnz()),
      g.work(res[0], this->nnz()), Tr, true, nrhs) << '\n';
  }

} // namespace casadi

#endif // CASADI_SOLVE_IMPL_HPP
