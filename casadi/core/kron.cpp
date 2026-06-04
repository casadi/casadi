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


#include "kron.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

namespace casadi {

  // ============ Sparsity-propagation traits (Fwd vs Rev) ============
  //
  // sp_forward and sp_reverse walk the exact same indices in the same order;
  // they differ only in the per-iteration bit-flow direction. Captured with
  // a `bool Fwd` template trait following the JacSparsityTraits pattern in
  // function_internal.cpp.

  template<bool Fwd> struct KronSpTraits;
  template<> struct KronSpTraits<true> {
    typedef const bvec_t** arg_t;
    typedef bvec_t** res_t;
    static inline void step(arg_t arg, res_t res,
                            casadi_int a_el, casadi_int b_el, casadi_int k) {
      res[0][k] = arg[0][a_el] | arg[1][b_el];
    }
  };
  template<> struct KronSpTraits<false> {
    typedef bvec_t** arg_t;
    typedef bvec_t** res_t;
    static inline void step(arg_t arg, res_t res,
                            casadi_int a_el, casadi_int b_el, casadi_int k) {
      arg[0][a_el] |= res[0][k];
      arg[1][b_el] |= res[0][k];
      res[0][k] = 0;
    }
  };

  // Templated Kron sparsity walk: identical 4-loop CSC scan in both modes.
  template<bool Fwd>
  static int kron_sp_gen(const Sparsity& sp_a, const Sparsity& sp_b,
                         typename KronSpTraits<Fwd>::arg_t arg,
                         typename KronSpTraits<Fwd>::res_t res) {
    const casadi_int* a_colind = sp_a.colind();
    const casadi_int* b_colind = sp_b.colind();
    const casadi_int a_ncol = sp_a.size2();
    const casadi_int b_ncol = sp_b.size2();
    casadi_int k = 0;
    for (casadi_int a_cc=0; a_cc<a_ncol; ++a_cc) {
      for (casadi_int b_cc=0; b_cc<b_ncol; ++b_cc) {
        for (casadi_int a_el=a_colind[a_cc]; a_el<a_colind[a_cc+1]; ++a_el) {
          for (casadi_int b_el=b_colind[b_cc]; b_el<b_colind[b_cc+1]; ++b_el) {
            KronSpTraits<Fwd>::step(arg, res, a_el, b_el, k);
            ++k;
          }
        }
      }
    }
    return 0;
  }


  // ============ Kron (base, both-sparse) ============

  MX Kron::create(const MX& a, const MX& b) {
    // Most-specific-first.
    if (auto* n = DenseKron::try_create(a, b))       return MX::create(n);
    if (auto* n = DenseSparseKron::try_create(a, b)) return MX::create(n);
    if (auto* n = SparseDenseKron::try_create(a, b)) return MX::create(n);
    return MX::create(new Kron(a, b));
  }

  Kron::Kron(const MX& a, const MX& b) {
    set_dep(a, b);
    set_sparsity(Sparsity::kron(a.sparsity(), b.sparsity()));
  }

  std::string Kron::disp(const std::vector<std::string>& arg) const {
    return "kron(" + arg.at(0) + ", " + arg.at(1) + ")";
  }

  void Kron::eval_kernel(const double** arg, double** res) const {
    casadi_kron(arg[0], dep(0).sparsity(), arg[1], dep(1).sparsity(), res[0]);
  }

  void Kron::eval_kernel(const SXElem** arg, SXElem** res) const {
    casadi_kron(arg[0], dep(0).sparsity(), arg[1], dep(1).sparsity(), res[0]);
  }

  void Kron::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const {
    res[0] = kron(arg[0], arg[1]);
  }

  int Kron::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return kron_sp_gen<true>(dep(0).sparsity(), dep(1).sparsity(), arg, res);
  }

  int Kron::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    return kron_sp_gen<false>(dep(0).sparsity(), dep(1).sparsity(), arg, res);
  }

  void Kron::ad_forward(const std::vector<std::vector<MX> >& fseed,
                        std::vector<std::vector<MX> >& fsens) const {
    const MX& A = dep(0);
    const MX& B = dep(1);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = kron(fseed[d][0], B) + kron(A, fseed[d][1]);
    }
  }

  void Kron::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                        std::vector<std::vector<MX> >& asens) const {
    const MX& A = dep(0);
    const MX& B = dep(1);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      const MX& Fbar = aseed[d][0];
      asens[d][0] += kron_contract(Fbar, B, true);
      asens[d][1] += kron_contract(Fbar, A, false);
    }
  }

  void Kron::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res,
                      const std::vector<bool>& arg_is_ref,
                      std::vector<bool>& res_is_ref) const {
    g.add_auxiliary(CodeGenerator::AUX_KRON);
    g << "casadi_kron("
      << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
      << g.sparsity(dep(0).sparsity()) << ", "
      << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
      << g.sparsity(dep(1).sparsity()) << ", "
      << g.work(res[0], nnz(), false) << ");\n";
  }

  void Kron::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Kron::kind", std::string("base"));
  }

  MXNode* Kron::deserialize(DeserializingStream& s) {
    std::string kind;
    s.unpack("Kron::kind", kind);
    if (kind == "base") return new Kron(s);
    if (kind == "dense") return new DenseKron(s);
    if (kind == "dense_sparse") return new DenseSparseKron(s);
    if (kind == "sparse_dense") return new SparseDenseKron(s);
    casadi_error("Unknown Kron kind: " + kind);
  }


  // ============ DenseKron ============

  MXNode* DenseKron::try_create(const MX& a, const MX& b) {
    if (!(a.is_dense() && b.is_dense())) return nullptr;
    return new DenseKron(a, b);
  }

  void DenseKron::eval_kernel(const double** arg, double** res) const {
    casadi_kron_dense(arg[0], dep(0).size1(), dep(0).size2(),
                      arg[1], dep(1).size1(), dep(1).size2(), res[0]);
  }

  void DenseKron::eval_kernel(const SXElem** arg, SXElem** res) const {
    casadi_kron_dense(arg[0], dep(0).size1(), dep(0).size2(),
                      arg[1], dep(1).size1(), dep(1).size2(), res[0]);
  }

  void DenseKron::generate(CodeGenerator& g,
                           const std::vector<casadi_int>& arg,
                           const std::vector<casadi_int>& res,
                           const std::vector<bool>& arg_is_ref,
                           std::vector<bool>& res_is_ref) const {
    g.add_auxiliary(CodeGenerator::AUX_KRON_DENSE);
    g << "casadi_kron_dense("
      << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
      << dep(0).size1() << ", " << dep(0).size2() << ", "
      << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
      << dep(1).size1() << ", " << dep(1).size2() << ", "
      << g.work(res[0], nnz(), false) << ");\n";
  }

  void DenseKron::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Kron::kind", std::string("dense"));
  }


  // ============ DenseSparseKron ============

  MXNode* DenseSparseKron::try_create(const MX& a, const MX& b) {
    if (!(a.is_dense() && !b.is_dense())) return nullptr;
    return new DenseSparseKron(a, b);
  }

  void DenseSparseKron::eval_kernel(const double** arg, double** res) const {
    casadi_kron_dense_sparse(arg[0], dep(0).size1(), dep(0).size2(),
                             arg[1], dep(1).sparsity(), res[0]);
  }

  void DenseSparseKron::eval_kernel(const SXElem** arg, SXElem** res) const {
    casadi_kron_dense_sparse(arg[0], dep(0).size1(), dep(0).size2(),
                             arg[1], dep(1).sparsity(), res[0]);
  }

  void DenseSparseKron::generate(CodeGenerator& g,
                                 const std::vector<casadi_int>& arg,
                                 const std::vector<casadi_int>& res,
                                 const std::vector<bool>& arg_is_ref,
                                 std::vector<bool>& res_is_ref) const {
    g.add_auxiliary(CodeGenerator::AUX_KRON_DENSE_SPARSE);
    g << "casadi_kron_dense_sparse("
      << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
      << dep(0).size1() << ", " << dep(0).size2() << ", "
      << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
      << g.sparsity(dep(1).sparsity()) << ", "
      << g.work(res[0], nnz(), false) << ");\n";
  }

  void DenseSparseKron::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Kron::kind", std::string("dense_sparse"));
  }


  // ============ SparseDenseKron ============

  MXNode* SparseDenseKron::try_create(const MX& a, const MX& b) {
    if (!(!a.is_dense() && b.is_dense())) return nullptr;
    return new SparseDenseKron(a, b);
  }

  void SparseDenseKron::eval_kernel(const double** arg, double** res) const {
    casadi_kron_sparse_dense(arg[0], dep(0).sparsity(),
                             arg[1], dep(1).size1(), dep(1).size2(), res[0]);
  }

  void SparseDenseKron::eval_kernel(const SXElem** arg, SXElem** res) const {
    casadi_kron_sparse_dense(arg[0], dep(0).sparsity(),
                             arg[1], dep(1).size1(), dep(1).size2(), res[0]);
  }

  void SparseDenseKron::generate(CodeGenerator& g,
                                 const std::vector<casadi_int>& arg,
                                 const std::vector<casadi_int>& res,
                                 const std::vector<bool>& arg_is_ref,
                                 std::vector<bool>& res_is_ref) const {
    g.add_auxiliary(CodeGenerator::AUX_KRON_SPARSE_DENSE);
    g << "casadi_kron_sparse_dense("
      << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
      << g.sparsity(dep(0).sparsity()) << ", "
      << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
      << dep(1).size1() << ", " << dep(1).size2() << ", "
      << g.work(res[0], nnz(), false) << ");\n";
  }

  void SparseDenseKron::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Kron::kind", std::string("sparse_dense"));
  }


  // ============ KronContract (base, fully general) ============

  MX KronContract::create(const MX& m, const MX& x, bool inner) {
    if (auto* n = DenseKronContract::try_create(m, x, inner))       return MX::create(n);
    if (auto* n = DenseSparseKronContract::try_create(m, x, inner)) return MX::create(n);
    if (auto* n = SparseDenseKronContract::try_create(m, x, inner)) return MX::create(n);
    return MX::create(new KronContract(m, x, inner));
  }

  KronContract::KronContract(const MX& m, const MX& x, bool inner) : inner_(inner) {
    casadi_assert(x.size1() > 0 && x.size2() > 0,
        "KronContract: X must be nonempty");
    casadi_assert(m.size1() % x.size1() == 0 && m.size2() % x.size2() == 0,
        "KronContract: M dims must be multiples of X dims");
    set_dep(m, x);
    set_sparsity(Sparsity::kron_contract(m.sparsity(), x.sparsity(), inner_));
  }

  std::string KronContract::disp(const std::vector<std::string>& arg) const {
    return std::string("kron_contract(") + arg.at(0) + ", " + arg.at(1)
         + ", " + (inner_ ? "inner" : "outer") + ")";
  }

  size_t KronContract::sz_w() const {
    return static_cast<size_t>(dep(1).numel());
  }

  void KronContract::eval_kernel(const double** arg, double** res, double* w) const {
    if (inner_) {
      casadi_kron_contract_inner(arg[0], dep(0).sparsity(),
                                 arg[1], dep(1).sparsity(),
                                 res[0], sparsity(), w);
    } else {
      casadi_kron_contract_outer(arg[0], dep(0).sparsity(),
                                 arg[1], dep(1).sparsity(),
                                 res[0], sparsity(), w);
    }
  }

  void KronContract::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    if (inner_) {
      casadi_kron_contract_inner(arg[0], dep(0).sparsity(),
                                 arg[1], dep(1).sparsity(),
                                 res[0], sparsity(), w);
    } else {
      casadi_kron_contract_outer(arg[0], dep(0).sparsity(),
                                 arg[1], dep(1).sparsity(),
                                 res[0], sparsity(), w);
    }
  }

  void KronContract::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const {
    res[0] = kron_contract(arg[0], arg[1], inner_);
  }

  int KronContract::sp_forward(const bvec_t** arg, bvec_t** res,
                               casadi_int* iw, bvec_t* w) const {
    // sp_forward and sp_reverse here have asymmetric X access (forward reads
    // X as input; reverse writes X-bar as output via per-iteration X-CSC
    // scan). Unifying via a shared template would require staging X-bar
    // through w in reverse and scattering at the end -- a real algorithmic
    // shift that can slow the (very sparse M, very dense X) corner by 100x+.
    // We keep them as two functions; the common scaffold is the column-
    // decomposition + row decomposition pattern which is short enough to
    // duplicate readably.
    const Sparsity& sp_m = dep(0).sparsity();
    const Sparsity& sp_x = dep(1).sparsity();
    const Sparsity& sp_y = sparsity();
    const casadi_int xrow = sp_x.size1(), xcol = sp_x.size2();
    const casadi_int yrow = sp_y.size1(), ycol = sp_y.size2();
    const casadi_int* m_colind = sp_m.colind();
    const casadi_int* m_row    = sp_m.row();
    const casadi_int* x_colind = sp_x.colind();
    const casadi_int* x_row    = sp_x.row();
    const casadi_int* y_colind = sp_y.colind();
    const casadi_int* y_row    = sp_y.row();
    const casadi_int mB = inner_ ? xrow : yrow;
    const casadi_int nB = inner_ ? xcol : ycol;
    const casadi_int nA = inner_ ? ycol : xcol;
    bvec_t* w_x = w;
    const casadi_int xn = xrow*xcol;
    for (casadi_int k=0; k<xn; ++k) w_x[k] = 0;
    for (casadi_int cc=0; cc<xcol; ++cc) {
      for (casadi_int el=x_colind[cc]; el<x_colind[cc+1]; ++el) {
        w_x[cc*xrow + x_row[el]] = arg[1][el];
      }
    }
    const casadi_int yn_total = y_colind[ycol];
    for (casadi_int k=0; k<yn_total; ++k) res[0][k] = 0;
    for (casadi_int j=0; j<nA; ++j) {
      for (casadi_int s=0; s<nB; ++s) {
        casadi_int cc = j*nB + s;
        casadi_int y_cc = inner_ ? j : s;
        casadi_int x_cc = inner_ ? s : j;
        casadi_int y_col_start = y_colind[y_cc];
        casadi_int y_col_end   = y_colind[y_cc+1];
        if (y_col_start == y_col_end) continue;
        for (casadi_int el=m_colind[cc]; el<m_colind[cc+1]; ++el) {
          casadi_int rr = m_row[el];
          casadi_int outer_row = rr / mB;
          casadi_int inner_row = rr % mB;
          casadi_int y_rr = inner_ ? outer_row : inner_row;
          casadi_int x_rr = inner_ ? inner_row : outer_row;
          bvec_t x_pat = w_x[x_cc*xrow + x_rr];
          for (casadi_int y_el=y_col_start; y_el<y_col_end; ++y_el) {
            if (y_row[y_el] == y_rr) { res[0][y_el] |= arg[0][el] | x_pat; break; }
            if (y_row[y_el] > y_rr) break;
          }
        }
      }
    }
    return 0;
  }

  int KronContract::sp_reverse(bvec_t** arg, bvec_t** res,
                               casadi_int* iw, bvec_t* w) const {
    const Sparsity& sp_m = dep(0).sparsity();
    const Sparsity& sp_x = dep(1).sparsity();
    const Sparsity& sp_y = sparsity();
    const casadi_int xrow = sp_x.size1(), xcol = sp_x.size2();
    const casadi_int yrow = sp_y.size1(), ycol = sp_y.size2();
    const casadi_int* m_colind = sp_m.colind();
    const casadi_int* m_row    = sp_m.row();
    const casadi_int* x_colind = sp_x.colind();
    const casadi_int* x_row    = sp_x.row();
    const casadi_int* y_colind = sp_y.colind();
    const casadi_int* y_row    = sp_y.row();
    const casadi_int mB = inner_ ? xrow : yrow;
    const casadi_int nB = inner_ ? xcol : ycol;
    const casadi_int nA = inner_ ? ycol : xcol;
    for (casadi_int j=0; j<nA; ++j) {
      for (casadi_int s=0; s<nB; ++s) {
        casadi_int cc = j*nB + s;
        casadi_int y_cc = inner_ ? j : s;
        casadi_int x_cc = inner_ ? s : j;
        casadi_int y_col_start = y_colind[y_cc];
        casadi_int y_col_end   = y_colind[y_cc+1];
        if (y_col_start == y_col_end) continue;
        casadi_int x_col_start = x_colind[x_cc];
        casadi_int x_col_end   = x_colind[x_cc+1];
        for (casadi_int el=m_colind[cc]; el<m_colind[cc+1]; ++el) {
          casadi_int rr = m_row[el];
          casadi_int outer_row = rr / mB;
          casadi_int inner_row = rr % mB;
          casadi_int y_rr = inner_ ? outer_row : inner_row;
          casadi_int x_rr = inner_ ? inner_row : outer_row;
          for (casadi_int y_el=y_col_start; y_el<y_col_end; ++y_el) {
            if (y_row[y_el] == y_rr) {
              bvec_t sval = res[0][y_el];
              arg[0][el] |= sval;
              for (casadi_int x_el=x_col_start; x_el<x_col_end; ++x_el) {
                if (x_row[x_el] == x_rr) { arg[1][x_el] |= sval; break; }
                if (x_row[x_el] > x_rr) break;
              }
              break;
            }
            if (y_row[y_el] > y_rr) break;
          }
        }
      }
    }
    const casadi_int yn_total = y_colind[ycol];
    for (casadi_int k=0; k<yn_total; ++k) res[0][k] = 0;
    return 0;
  }

  void KronContract::ad_forward(const std::vector<std::vector<MX> >& fseed,
                                std::vector<std::vector<MX> >& fsens) const {
    const MX& M = dep(0);
    const MX& X = dep(1);
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = kron_contract(fseed[d][0], X, inner_)
                  + kron_contract(M, fseed[d][1], inner_);
    }
  }

  void KronContract::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                                std::vector<std::vector<MX> >& asens) const {
    const MX& M = dep(0);
    const MX& X = dep(1);
    for (casadi_int d=0; d<aseed.size(); ++d) {
      const MX& Ybar = aseed[d][0];
      MX m_kron = inner_ ? kron(Ybar, X) : kron(X, Ybar);
      asens[d][0] += project(m_kron, M.sparsity());
      asens[d][1] += kron_contract(M, Ybar, !inner_);
    }
  }

  void KronContract::generate(CodeGenerator& g,
                              const std::vector<casadi_int>& arg,
                              const std::vector<casadi_int>& res,
                              const std::vector<bool>& arg_is_ref,
                              std::vector<bool>& res_is_ref) const {
    if (inner_) {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_INNER);
      g << "casadi_kron_contract_inner(";
    } else {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_OUTER);
      g << "casadi_kron_contract_outer(";
    }
    g << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
      << g.sparsity(dep(0).sparsity()) << ", "
      << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
      << g.sparsity(dep(1).sparsity()) << ", "
      << g.work(res[0], nnz(), false) << ", "
      << g.sparsity(sparsity()) << ", w);\n";
  }

  void KronContract::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("KronContract::inner", inner_);
  }

  void KronContract::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("KronContract::kind", std::string("base"));
  }

  KronContract::KronContract(DeserializingStream& s) : MXNode(s) {
    s.unpack("KronContract::inner", inner_);
  }

  MXNode* KronContract::deserialize(DeserializingStream& s) {
    std::string kind;
    s.unpack("KronContract::kind", kind);
    if (kind == "base") return new KronContract(s);
    if (kind == "dense") return new DenseKronContract(s);
    if (kind == "dense_sparse") return new DenseSparseKronContract(s);
    if (kind == "sparse_dense") return new SparseDenseKronContract(s);
    casadi_error("Unknown KronContract kind: " + kind);
  }


  // ============ DenseKronContract ============

  MXNode* DenseKronContract::try_create(const MX& m, const MX& x, bool inner) {
    if (!(m.is_dense() && x.is_dense())) return nullptr;
    return new DenseKronContract(m, x, inner);
  }

  void DenseKronContract::eval_kernel(const double** arg, double** res, double* w) const {
    if (inner_) {
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      const casadi_int mB = dep(1).size1(),     nB = dep(1).size2();
      casadi_kron_contract_inner_dense(arg[0], mA, nA, arg[1], mB, nB, res[0]);
    } else {
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      const casadi_int mA = dep(1).size1(),     nA = dep(1).size2();
      casadi_kron_contract_outer_dense(arg[0], mB, nB, arg[1], mA, nA, res[0]);
    }
  }

  void DenseKronContract::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    if (inner_) {
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      const casadi_int mB = dep(1).size1(),     nB = dep(1).size2();
      casadi_kron_contract_inner_dense(arg[0], mA, nA, arg[1], mB, nB, res[0]);
    } else {
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      const casadi_int mA = dep(1).size1(),     nA = dep(1).size2();
      casadi_kron_contract_outer_dense(arg[0], mB, nB, arg[1], mA, nA, res[0]);
    }
  }

  void DenseKronContract::generate(CodeGenerator& g,
                                   const std::vector<casadi_int>& arg,
                                   const std::vector<casadi_int>& res,
                                   const std::vector<bool>& arg_is_ref,
                                   std::vector<bool>& res_is_ref) const {
    if (inner_) {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_INNER_DENSE);
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      const casadi_int mB = dep(1).size1(),     nB = dep(1).size2();
      g << "casadi_kron_contract_inner_dense("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << mA << ", " << nA << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << mB << ", " << nB << ", "
        << g.work(res[0], nnz(), false) << ");\n";
    } else {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_OUTER_DENSE);
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      const casadi_int mA = dep(1).size1(),     nA = dep(1).size2();
      g << "casadi_kron_contract_outer_dense("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << mB << ", " << nB << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << mA << ", " << nA << ", "
        << g.work(res[0], nnz(), false) << ");\n";
    }
  }

  void DenseKronContract::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("KronContract::kind", std::string("dense"));
  }


  // ============ DenseSparseKronContract ============

  MXNode* DenseSparseKronContract::try_create(const MX& m, const MX& x, bool inner) {
    if (!(m.is_dense() && !x.is_dense())) return nullptr;
    return new DenseSparseKronContract(m, x, inner);
  }

  void DenseSparseKronContract::eval_kernel(const double** arg, double** res, double* w) const {
    if (inner_) {
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      casadi_kron_contract_inner_dense_sparse(arg[0], mA, nA,
                                              arg[1], dep(1).sparsity(), res[0]);
    } else {
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      casadi_kron_contract_outer_dense_sparse(arg[0], mB, nB,
                                              arg[1], dep(1).sparsity(), res[0]);
    }
  }

  void DenseSparseKronContract::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    if (inner_) {
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      casadi_kron_contract_inner_dense_sparse(arg[0], mA, nA,
                                              arg[1], dep(1).sparsity(), res[0]);
    } else {
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      casadi_kron_contract_outer_dense_sparse(arg[0], mB, nB,
                                              arg[1], dep(1).sparsity(), res[0]);
    }
  }

  void DenseSparseKronContract::generate(CodeGenerator& g,
                                         const std::vector<casadi_int>& arg,
                                         const std::vector<casadi_int>& res,
                                         const std::vector<bool>& arg_is_ref,
                                         std::vector<bool>& res_is_ref) const {
    if (inner_) {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_INNER_DENSE_SPARSE);
      const casadi_int mA = sparsity().size1(), nA = sparsity().size2();
      g << "casadi_kron_contract_inner_dense_sparse("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << mA << ", " << nA << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << g.sparsity(dep(1).sparsity()) << ", "
        << g.work(res[0], nnz(), false) << ");\n";
    } else {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_OUTER_DENSE_SPARSE);
      const casadi_int mB = sparsity().size1(), nB = sparsity().size2();
      g << "casadi_kron_contract_outer_dense_sparse("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << mB << ", " << nB << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << g.sparsity(dep(1).sparsity()) << ", "
        << g.work(res[0], nnz(), false) << ");\n";
    }
  }

  void DenseSparseKronContract::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("KronContract::kind", std::string("dense_sparse"));
  }


  // ============ SparseDenseKronContract ============

  MXNode* SparseDenseKronContract::try_create(const MX& m, const MX& x, bool inner) {
    if (!(!m.is_dense() && x.is_dense())) return nullptr;
    return new SparseDenseKronContract(m, x, inner);
  }

  void SparseDenseKronContract::eval_kernel(const double** arg, double** res, double* w) const {
    if (inner_) {
      casadi_kron_contract_inner_sparse_dense(arg[0], dep(0).sparsity(),
                                              arg[1], dep(1).size1(), dep(1).size2(),
                                              res[0], sparsity());
    } else {
      casadi_kron_contract_outer_sparse_dense(arg[0], dep(0).sparsity(),
                                              arg[1], dep(1).size1(), dep(1).size2(),
                                              res[0], sparsity());
    }
  }

  void SparseDenseKronContract::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    if (inner_) {
      casadi_kron_contract_inner_sparse_dense(arg[0], dep(0).sparsity(),
                                              arg[1], dep(1).size1(), dep(1).size2(),
                                              res[0], sparsity());
    } else {
      casadi_kron_contract_outer_sparse_dense(arg[0], dep(0).sparsity(),
                                              arg[1], dep(1).size1(), dep(1).size2(),
                                              res[0], sparsity());
    }
  }

  void SparseDenseKronContract::generate(CodeGenerator& g,
                                         const std::vector<casadi_int>& arg,
                                         const std::vector<casadi_int>& res,
                                         const std::vector<bool>& arg_is_ref,
                                         std::vector<bool>& res_is_ref) const {
    if (inner_) {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_INNER_SPARSE_DENSE);
      g << "casadi_kron_contract_inner_sparse_dense("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << g.sparsity(dep(0).sparsity()) << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << dep(1).size1() << ", " << dep(1).size2() << ", "
        << g.work(res[0], nnz(), false) << ", "
        << g.sparsity(sparsity()) << ");\n";
    } else {
      g.add_auxiliary(CodeGenerator::AUX_KRON_CONTRACT_OUTER_SPARSE_DENSE);
      g << "casadi_kron_contract_outer_sparse_dense("
        << g.work(arg[0], dep(0).nnz(), arg_is_ref[0]) << ", "
        << g.sparsity(dep(0).sparsity()) << ", "
        << g.work(arg[1], dep(1).nnz(), arg_is_ref[1]) << ", "
        << dep(1).size1() << ", " << dep(1).size2() << ", "
        << g.work(res[0], nnz(), false) << ", "
        << g.sparsity(sparsity()) << ");\n";
    }
  }

  void SparseDenseKronContract::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("KronContract::kind", std::string("sparse_dense"));
  }

} // namespace casadi
