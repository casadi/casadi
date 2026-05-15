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


#ifndef CASADI_MULTIPLICATION_CPP
#define CASADI_MULTIPLICATION_CPP

#include "multiplication.hpp"
#include "blas_impl.hpp"
#include "casadi_misc.hpp"
#include "function_internal.hpp"
#include "serializing_stream.hpp"

namespace casadi {

  MX Multiplication::create(const MX& z, const MX& x, const MX& y,
                            const std::string& blas) {
    // Most-specific-first. Each subclass owns its own applicability test
    // and construction; nullptr means "doesn't apply, try the next one".
    if (auto* n = DenseMultiplication::try_create(z, x, y, blas))       return MX::create(n);
    if (auto* n = PseudoDenseMultiplication::try_create(z, x, y, blas)) return MX::create(n);
    if (auto* n = DenseSparseMultiplication::try_create(z, x, y, blas)) return MX::create(n);
    return MX::create(new Multiplication(z, x, y, blas));
  }

  Multiplication::Multiplication(const MX& z, const MX& x, const MX& y,
                                 const std::string& blas) {
    casadi_assert(x.size2() == y.size1() && x.size1() == z.size1()
      && y.size2() == z.size2(),
      "Multiplication::Multiplication: dimension mismatch. Attempting to multiply "
      + x.dim() + " with " + y.dim()
      + " and add the result to " + z.dim());

    set_dep(z, x, y);
    set_sparsity(z.sparsity());
    blas_shorthand_ = Blas::shorthand_for(blas);
  }

  std::string Multiplication::disp(const std::vector<std::string>& arg) const {
    return "mac(" + arg.at(1) + "," + arg.at(2) + "," + arg.at(0) + ")";
  }

  void Multiplication::eval_kernel(const double** arg, double** res, double* w) const {
    casadi_mtimes(arg[1], dep(1).sparsity(), arg[2], dep(2).sparsity(),
                  res[0], sparsity(), w, false);
  }

  void Multiplication::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    casadi_mtimes(arg[1], dep(1).sparsity(), arg[2], dep(2).sparsity(),
                  res[0], sparsity(), w, false);
  }

  void Multiplication::ad_forward(const std::vector<std::vector<MX> >& fseed,
                               std::vector<std::vector<MX> >& fsens) const {
    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = fseed[d][0]
        + mac(dep(1), fseed[d][2], MX::zeros(dep(0).sparsity()))
        + mac(fseed[d][1], dep(2), MX::zeros(dep(0).sparsity()));
    }
  }

  void Multiplication::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                               std::vector<std::vector<MX> >& asens) const {
    for (casadi_int d=0; d<aseed.size(); ++d) {
      MX dep2T = dep(2).is_constant() ? MX(DM(dep(2)).T()) : dep(2).T();
      asens[d][1] += mac(aseed[d][0], dep2T, MX::zeros(dep(1).sparsity()));
      MX dep1T = dep(1).is_constant() ? MX(DM(dep(1)).T()) : dep(1).T();
      asens[d][2] += mac(dep1T, aseed[d][0], MX::zeros(dep(2).sparsity()));
      asens[d][0] += aseed[d][0];
    }
  }

  void Multiplication::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const {
    res[0] = mac(arg[1], arg[2], arg[0]);
  }

  void Multiplication::eval_linear(const std::vector<std::array<MX, 3> >& arg,
                        std::vector<std::array<MX, 3> >& res) const {
    const std::array<MX, 3>& x = arg[1];
    const std::array<MX, 3>& y = arg[2];
    const std::array<MX, 3>& z = arg[0];
    std::array<MX, 3>& f = res[0];
    f[0] = mac(x[0], y[0], z[0]);
    f[1] = mac(x[0], y[1], z[1]);
    f[1] = mac(x[1], y[0], f[1]);
    f[2] = mac(x[1]+x[2], y[1]+y[2], z[2]);
    f[2] = mac(x[2], y[0], f[2]);
    f[2] = mac(x[0], y[2], f[2]);
  }

  int Multiplication::
  sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    copy_fwd(arg[0], res[0], nnz());
    Sparsity::mul_sparsityF(arg[1], dep(1).sparsity(),
                            arg[2], dep(2).sparsity(),
                            res[0], sparsity(), w);
    return 0;
  }

  int Multiplication::
  sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    Sparsity::mul_sparsityR(arg[1], dep(1).sparsity(),
                            arg[2], dep(2).sparsity(),
                            res[0], sparsity(), w);
    copy_rev(arg[0], res[0], nnz());
    return 0;
  }

  // Helper: emit `if (arg!=res || arg_is_ref) copy(arg -> res, nnz)`
  // Subclasses share this so eval_kernel only handles z += x*y.
  static void codegen_copy_z(CodeGenerator& g, casadi_int nnz,
                             const std::vector<casadi_int>& arg,
                             const std::vector<casadi_int>& res,
                             const std::vector<bool>& arg_is_ref) {
    if (arg[0]!=res[0] || arg_is_ref[0]) {
      g << g.copy(g.work(arg[0], nnz, arg_is_ref[0]),
                  nnz,
                  g.work(res[0], nnz, false)) << '\n';
    }
  }

  void Multiplication::generate(CodeGenerator& g,
                                const std::vector<casadi_int>& arg,
                                const std::vector<casadi_int>& res,
                                const std::vector<bool>& arg_is_ref,
                                std::vector<bool>& res_is_ref) const {
    codegen_copy_z(g, nnz(), arg, res, arg_is_ref);
    g << g.mtimes(g.work(arg[1], dep(1).nnz(), arg_is_ref[1]), dep(1).sparsity(),
                  g.work(arg[2], dep(2).nnz(), arg_is_ref[2]), dep(2).sparsity(),
                  g.work(res[0], nnz(), false), sparsity(), "w", false) << '\n';
  }

  // ---------------- DenseMultiplication ----------------

  MXNode* DenseMultiplication::try_create(const MX& z, const MX& x, const MX& y,
                                          const std::string& blas) {
    if (!(x.is_dense() && y.is_dense() && z.is_dense())) return nullptr;
    return new DenseMultiplication(z, x, y, blas);
  }

  void DenseMultiplication::eval_kernel(const double** arg, double** res, double* w) const {
    Blas::mtimes(blas_shorthand_,
                 arg[1], dep(1).size1(), dep(1).size2(),
                 arg[2], dep(2).size2(),
                 res[0]);
  }

  void DenseMultiplication::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    casadi_mtimes_dense(arg[1], dep(1).size1(), dep(1).size2(),
                        arg[2], dep(2).size2(), res[0], false);
  }

  void DenseMultiplication::generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res,
           const std::vector<bool>& arg_is_ref, std::vector<bool>& res_is_ref) const {
    codegen_copy_z(g, nnz(), arg, res, arg_is_ref);
    Blas::codegen_mtimes(g, blas_shorthand_,
        g.work(arg[1], dep(1).nnz(), arg_is_ref[1]),
        dep(1).size1(), dep(1).size2(),
        g.work(arg[2], dep(2).nnz(), arg_is_ref[2]),
        dep(2).size2(),
        g.work(res[0], nnz(), false));
  }

  // ---------------- DenseSparseMultiplication ----------------

  MXNode* DenseSparseMultiplication::try_create(const MX& z, const MX& x, const MX& y,
                                                const std::string& blas) {
    if (!(x.is_dense() && z.is_dense())) return nullptr;
    return new DenseSparseMultiplication(z, x, y, blas);
  }

  void DenseSparseMultiplication::eval_kernel(const double** arg, double** res, double* w) const {
    casadi_mtimes_dense_sparse(arg[1], dep(1).size1(),
                               arg[2], dep(2).sparsity(), res[0]);
  }

  void DenseSparseMultiplication::eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const {
    casadi_mtimes_dense_sparse(arg[1], dep(1).size1(),
                               arg[2], dep(2).sparsity(), res[0]);
  }

  void DenseSparseMultiplication::generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res,
           const std::vector<bool>& arg_is_ref, std::vector<bool>& res_is_ref) const {
    codegen_copy_z(g, nnz(), arg, res, arg_is_ref);
    g << g.mtimes_dense_sparse(g.work(arg[1], dep(1).nnz(), arg_is_ref[1]), dep(1).size1(),
                               g.work(arg[2], dep(2).nnz(), arg_is_ref[2]), dep(2).sparsity(),
                               g.work(res[0], nnz(), false)) << '\n';
  }

  // ---------------- PseudoDenseMultiplication ----------------

  MXNode* PseudoDenseMultiplication::try_create(const MX& z, const MX& x, const MX& y,
                                                const std::string& blas) {
    // Result must have at least one nonzero for the compact-buffer reinterpretation
    // to be meaningful.
    if (z.nnz() == 0) return nullptr;
    std::vector<casadi_int> xr, xc, yr, yc, zr, zc;
    if (!x.sparsity().is_compactible(xr, xc)) return nullptr;
    if (!y.sparsity().is_compactible(yr, yc)) return nullptr;
    if (!z.sparsity().is_compactible(zr, zc)) return nullptr;
    // Connecting index sets must agree exactly: cols(x) = rows(y), and
    // result rows/cols inherit from x/y.
    if (xc != yr || xr != zr || yc != zc) return nullptr;
    return new PseudoDenseMultiplication(z, x, y,
        static_cast<casadi_int>(xr.size()),
        static_cast<casadi_int>(xc.size()),
        static_cast<casadi_int>(yc.size()),
        blas);
  }

  PseudoDenseMultiplication::PseudoDenseMultiplication(const MX& z, const MX& x, const MX& y,
      casadi_int nrow_x_compact, casadi_int ncol_x_compact, casadi_int ncol_y_compact,
      const std::string& blas)
    : Multiplication(z, x, y, blas),
      a_(nrow_x_compact), b_(ncol_x_compact), c_(ncol_y_compact) {}

  void PseudoDenseMultiplication::eval_kernel(const double** arg,
      double** res, double* w) const {
    Blas::mtimes(blas_shorthand_, arg[1], a_, b_, arg[2], c_, res[0]);
  }

  void PseudoDenseMultiplication::eval_kernel(const SXElem** arg,
      SXElem** res, SXElem* w) const {
    casadi_mtimes_dense(arg[1], a_, b_, arg[2], c_, res[0], false);
  }

  void PseudoDenseMultiplication::generate(CodeGenerator& g,
           const std::vector<casadi_int>& arg, const std::vector<casadi_int>& res,
           const std::vector<bool>& arg_is_ref, std::vector<bool>& res_is_ref) const {
    codegen_copy_z(g, nnz(), arg, res, arg_is_ref);
    Blas::codegen_mtimes(g, blas_shorthand_,
        g.work(arg[1], dep(1).nnz(), arg_is_ref[1]), a_, b_,
        g.work(arg[2], dep(2).nnz(), arg_is_ref[2]), c_,
        g.work(res[0], nnz(), false));
  }

  void Multiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Multiplication::kind", std::string("base"));
  }

  void Multiplication::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("Multiplication::blas",
           std::string(Blas::name_for_shorthand(blas_shorthand_)));
  }

  Multiplication::Multiplication(DeserializingStream& s, bool legacy) : MXNode(s) {
    if (legacy) {
      blas_shorthand_ = 0;
      return;
    }
    std::string blas;
    s.unpack("Multiplication::blas", blas);
    blas_shorthand_ = Blas::shorthand_for(blas);
  }

  void DenseMultiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Multiplication::kind", std::string("dense"));
  }

  void DenseSparseMultiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Multiplication::kind", std::string("dense_sparse"));
  }

  void PseudoDenseMultiplication::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("Multiplication::kind", std::string("pseudo_dense"));
  }

  void PseudoDenseMultiplication::serialize_body(SerializingStream& s) const {
    Multiplication::serialize_body(s);
    s.pack("PseudoDenseMultiplication::a", a_);
    s.pack("PseudoDenseMultiplication::b", b_);
    s.pack("PseudoDenseMultiplication::c", c_);
  }

  PseudoDenseMultiplication::PseudoDenseMultiplication(DeserializingStream& s)
    : Multiplication(s) {
    s.unpack("PseudoDenseMultiplication::a", a_);
    s.unpack("PseudoDenseMultiplication::b", b_);
    s.unpack("PseudoDenseMultiplication::c", c_);
  }

  MXNode* Multiplication::deserialize(DeserializingStream& s) {
    // Wire-format detection.
    //   pre 3.8: serialize_type packed a *bool* ("Multiplication::dense").
    //     Body did NOT include a blas field. Two variants only --
    //     DenseMultiplication and the generic Multiplication.
    //   3.8+:    serialize_type packs a *string* ("Multiplication::kind"),
    //     four variants, body always includes "Multiplication::blas".
    //
    // In debug mode the descriptor is on the wire and IS the discriminator;
    // in non-debug mode we discriminate on the first wire byte: bool encodes
    // as 2 hex chars whose first is 'a' or 'b'; string-length first byte is
    // 'a'+(length%16), which for our kind names {"base","dense",
    // "dense_sparse","pseudo_dense"} (lengths 4,5,12,12) is 'e','f','m','m'.
    bool legacy_bool;
    std::string kind;
    if (s.debug()) {
      std::string descr;
      s.unpack(descr);
      if (descr == "Multiplication::dense") {
        legacy_bool = true;
      } else {
        casadi_assert(descr == "Multiplication::kind",
            "Unexpected Multiplication descriptor: '" + descr + "'.");
        legacy_bool = false;
        s.unpack(kind);
      }
    } else {
      const int p = s.peek_byte();
      legacy_bool = (p == 'a' || p == 'b');
      if (!legacy_bool) s.unpack("Multiplication::kind", kind);
    }

    if (legacy_bool) {
      bool dense;
      if (s.debug()) s.unpack(dense);
      else           s.unpack("Multiplication::dense", dense);
      if (dense) return new DenseMultiplication(s, /*legacy=*/true);
      return new Multiplication(s, /*legacy=*/true);
    }

    if (kind == "dense")        return new DenseMultiplication(s);
    if (kind == "dense_sparse") return new DenseSparseMultiplication(s);
    if (kind == "pseudo_dense") return new PseudoDenseMultiplication(s);
    return new Multiplication(s);
  }

} // namespace casadi

#endif // CASADI_MULTIPLICATION_CPP
