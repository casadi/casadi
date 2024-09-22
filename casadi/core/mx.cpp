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


#include "mx_node.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "multiple_output.hpp"
#include "casadi_misc.hpp"
#include "norm.hpp"
#include "calculus.hpp"
#include "mx_function.hpp"
#include "linsol.hpp"
#include "expm.hpp"
#include "serializing_stream.hpp"
#include "im.hpp"
#include "bspline.hpp"
#include "casadi_call.hpp"

// Throw informative error message
#define CASADI_THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

// Throw informative error message
#define CASADI_THROW_ERROR_OBJ(FNAME, WHAT) \
throw CasadiException("Error in MX::" FNAME " for node of type " \
  + this->class_name() + " at " + CASADI_WHERE + ":\n" + std::string(WHAT));

namespace casadi {

  template class GenericMatrix< MX >;

  MX::~MX() {
  }

  MX::MX() {
    own(ZeroByZero::getInstance());
  }

  MX::MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4) {
    own(node);
  }

  MX MX::create(MXNode* node) {
    return MX(node, false, false, false, false);
  }

  MX::MX(double x) {
    own(ConstantMX::create(Sparsity::dense(1, 1), x));
  }

  MX::MX(const DM& x) {
    own(ConstantMX::create(x));
  }

  MX::MX(const std::vector<double>& x) {
    own(ConstantMX::create(DM(x)));
  }

  MX::MX(const Sparsity& sp, const MX& val) {
    if (sp.is_reshape(val.sparsity())) {
      *this = reshape(val, sp);
    } else if (val.is_scalar()) {
      // Dense matrix if val dense
      if (val.is_dense()) {
        if (val.is_constant()) {
          own(ConstantMX::create(sp, static_cast<double>(val)));
        } else {
          *this = val->get_nzref(sp, std::vector<casadi_int>(sp.nnz(), 0));
        }
      } else {
        // Empty matrix
        own(ConstantMX::create(Sparsity(sp.size()), 0));
      }
    } else {
      casadi_assert_dev(val.is_column() && sp.nnz()==val.size1());
      *this = densify(val)->get_nzref(sp, range(sp.nnz()));
    }
  }

  MX::MX(const Sparsity& sp) {
    own(ConstantMX::create(sp, 1));
  }

  MX::MX(casadi_int nrow, casadi_int ncol) {
    own(ConstantMX::create(Sparsity(nrow, ncol), 0));
  }

  MX::MX(const std::pair<casadi_int, casadi_int>& rc) {
    own(ConstantMX::create(Sparsity(rc), 0));
  }

  MX::MX(const Sparsity& sp, double val, bool dummy) {
    own(ConstantMX::create(sp, val));
  }

  MX::MX(const Sparsity& sp, const std::string& fname) {
    own(ConstantMX::create(sp, fname));
  }

  MX::MX(const DM& val, const std::string& name) {
    own(ConstantMX::create(val, name));
  }

  std::vector<MX> MX::createMultipleOutput(MXNode* node) {
    casadi_assert_dev(dynamic_cast<MultipleOutput*>(node) != nullptr);
    MX x =  MX::create(node);
    std::vector<MX> ret(x->nout());
    for (casadi_int i=0; i<ret.size(); ++i) {
      ret[i] = x.get_output(i);
      if (ret[i].is_empty(true)) {
        ret[i] = MX(0, 0);
      } else if (ret[i].nnz()==0) {
        ret[i] = MX(ret[i].size());
      }
    }
    return ret;
  }

  bool MX::__nonzero__() const {
    return (*this)->__nonzero__();
  }

  void MX::get(MX& m, bool ind1, const Slice& rr, const Slice& cc) const {
    // Fall back on (IM, IM)
    return get(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc) const {
    // Fall back on (IM, IM)
    get(m, ind1, rr.all(size1(), ind1), cc);
  }

  void MX::get(MX& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc) const {
    // Fall back on (IM, IM)
    get(m, ind1, rr, cc.all(size2(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) const {
    // Make sure dense vectors
    casadi_assert(rr.is_dense() && rr.is_vector(),
                          "Marix::get: First index must be a dense vector");
    casadi_assert(cc.is_dense() && cc.is_vector(),
                          "Marix::get: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), cc.nonzeros(), mapping, ind1);

    // Create return MX
    m = (*this)->get_nzref(sp, mapping);
  }

  void MX::get(MX& m, bool ind1, const Slice& rr) const {
    // Fall back on IM
    get(m, ind1, rr.all(numel(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Matrix<casadi_int>& rr) const {
    // If the indexed matrix is dense, use nonzero indexing
    if (is_dense()) {
      return get_nz(m, ind1, rr);
    }

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && rr.is_row()) || (is_row() && rr.is_column());

    // Get the sparsity pattern - does bounds checking
    std::vector<casadi_int> mapping;
    Sparsity sp = sparsity().sub(rr.nonzeros(), tr ? rr.sparsity().T() : rr.sparsity(),
                                 mapping, ind1);

    // Create return MX
    m = (*this)->get_nzref(sp, mapping);
  }

  void MX::get(MX& m, bool ind1, const Sparsity& sp) const {
    casadi_assert(size()==sp.size(),
      "get(Sparsity sp): shape mismatch. This matrix has shape "
      + str(size()) + ", but supplied sparsity index has shape "
      + str(sp.size()) + ".");
    m = project(*this, sp);
  }

  void MX::get(MX& m, bool ind1, const MX& rr) const {
    casadi_assert(is_dense(), "Parametric slicing only supported for dense matrices."
                              "Got " + dim(true) + " instead.");
    get_nz(m, ind1, rr);
  }

  void MX::get(MX& m, bool ind1, const Slice& rr, const MX& cc) const {
    casadi_assert(is_dense(), "Parametric slicing only supported for dense matrices. ");
    m = (*this)->get_nz_ref(rr.apply(size1()), floor(ind1 ? cc-1 : cc)*size1());
  }

  void MX::get(MX& m, bool ind1, const MX& rr, const Slice& cc) const {
    casadi_assert(is_dense(), "Parametric slicing only supported for dense matrices.");
    m = (*this)->get_nz_ref(ind1 ? rr-1 : rr, cc.apply(size2())*size1());
  }

  void MX::get(MX& m, bool ind1, const MX& rr, const MX& cc) const {
    casadi_assert(is_dense(), "Parametric slicing only supported for dense matrices.");
    m = (*this)->get_nz_ref(ind1 ? rr-1 : rr, floor(ind1 ? cc-1 : cc)*size1());
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc.all(size2(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr, const Matrix<casadi_int>& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.all(size1(), ind1), cc);
  }

  void MX::set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr, cc.all(size2(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Matrix<casadi_int>& rr, const Matrix<casadi_int>& cc) {
    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return set(m, ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return set(m, ind1, rr, cc.T());
    }

    // Make sure rr and cc are dense vectors
    casadi_assert(rr.is_dense() && rr.is_column(),
      "MX::set: First index not dense vector");
    casadi_assert(cc.is_dense() && cc.is_column(),
      "MX::set: Second index not dense vector");

    // Assert dimensions of assigning matrix
    if (rr.size1() != m.size1() || cc.size1() != m.size2()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        return set(repmat(m, rr.size1(), cc.size1()), ind1, rr, cc);
      } else if (rr.size1() == m.size2() && cc.size1() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr, cc);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(rr.size1()) + "-by-"
                     + str(cc.size1()) + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions
    casadi_int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    casadi_assert_in_range(rr.nonzeros(), -sz1+ind1, sz1+ind1);
    casadi_assert_in_range(cc.nonzeros(), -sz2+ind1, sz2+ind1);

    // If we are assigning with something sparse, first remove existing entries
    if (!m.is_dense()) {
      erase(rr.nonzeros(), cc.nonzeros(), ind1);
    }

    // Collect all assignments
    IM el = IM::zeros(m.sparsity());
    for (casadi_int j=0; j<el.size2(); ++j) { // Loop over columns of m
      casadi_int this_j = cc->at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (casadi_int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        casadi_int i = m.row(k);
        casadi_int this_i = rr->at(i) - ind1; // Corresponding row in this
        if (this_i<0) this_i += sz1;
        el->at(k) = this_i + this_j*sz1;
      }
    }
    return set(m, false, el);
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr) {
    // Fall back on IM
    set(m, ind1, rr.all(size1(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Matrix<casadi_int>& rr) {
    // Assert dimensions of assigning matrix
    if (rr.sparsity() != m.sparsity()) {
      if (rr.size() == m.size()) {
        // Remove submatrix to be replaced
        erase(rr.nonzeros(), ind1);

        // Find the intersection between rr's and m's sparsity patterns
        Sparsity sp = rr.sparsity() * m.sparsity();

        // Project both matrices to this sparsity
        return set(project(m, sp), ind1, Matrix<casadi_int>::project(rr, sp));
      } else if (m.is_scalar()) {
        // m scalar means "set all"
        if (m.is_dense()) {
          return set(MX(rr.sparsity(), m), ind1, rr);
        } else {
          return set(MX(rr.size()), ind1, rr);
        }
      } else if (rr.size1() == m.size2() && rr.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set(m.T(), ind1, rr);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(rr.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Dimensions of this
    casadi_int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    casadi_assert_in_range(rr.nonzeros(), -nel+ind1, nel+ind1);

    // Dense mode
    if (is_dense() && m.is_dense()) {
      return set_nz(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<casadi_int> new_row=sparsity().get_row(), new_col=sparsity().get_col();
    std::vector<casadi_int> nz(rr.nonzeros());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    for (std::vector<casadi_int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (ind1) (*i)--;
      if (*i<0) *i += nel;
      new_row.push_back(*i % sz1);
      new_col.push_back(*i / sz1);
    }
    Sparsity sp = Sparsity::triplet(sz1, sz2, new_row, new_col);

    // If needed, update pattern
    if (sp != sparsity()) *this = project(*this, sp);

    // Find the nonzeros corresponding to rr
    sparsity().get_nz(nz);

    // Create a nonzero assignment node
    *this = m->get_nzassign(*this, nz);
  }

  void MX::set(const MX& m, bool ind1, const Sparsity& sp) {
    casadi_assert(size()==sp.size(),
      "set(Sparsity sp): shape mismatch. This matrix has shape "
      + str(size()) + ", but supplied sparsity index has shape "
      + str(sp.size()) + ".");
    std::vector<casadi_int> ii = sp.find();
    if (m.is_scalar()) {
      (*this)(ii) = densify(m);
    } else {
      (*this)(ii) = densify(m(ii));
    }
  }

  void MX::get_nz(MX& m, bool ind1, const Slice& kk) const {
    // Fallback on IM
    get_nz(m, ind1, kk.all(nnz(), ind1));
  }

  void MX::get_nz(MX& m, bool ind1, const Matrix<casadi_int>& kk) const {
    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && kk.is_row()) || (is_row() && kk.is_column());

    // Quick return if no entries
    if (kk.nnz()==0) {
      m = MX::zeros(tr ? kk.sparsity().T() : kk.sparsity());
      return;
    }

    // Check bounds
    casadi_int sz = nnz();
    casadi_assert_in_range(kk.nonzeros(), -sz+ind1, sz+ind1);

    // Handle index-1, negative indices
    if (ind1 || *std::min_element(kk->begin(), kk->end())<0) {
      Matrix<casadi_int> kk_mod = kk;
      for (auto&& i : kk_mod.nonzeros()) {
        casadi_assert(!(ind1 && i<=0),
          "Matlab is 1-based, but requested index " + str(i) + ". "
          "Note that negative slices are disabled in the Matlab interface. "
          "Possibly you may want to use 'end'.");
        if (ind1) i--;
        if (i<0) i += sz;
      }
      get_nz(m, false, kk_mod); // Call recursively
      return;
    }

    // Return reference to the nonzeros
    m = (*this)->get_nzref(tr ? kk.sparsity().T() : kk.sparsity(), kk.nonzeros());
  }

  void MX::get_nz(MX& m, bool ind1, const MX& kk) const {
    // Create return MX
    m = (*this)->get_nz_ref(ind1 ? kk-1.0 : kk);
  }

  void MX::get_nz(MX& m, bool ind1, const MX& inner, const MX& outer) const {
    // Create return MX
    m = (*this)->get_nz_ref(ind1 ? inner-1.0: inner, ind1 ? outer-1.0: outer);
  }

  void MX::get_nz(MX& m, bool ind1, const Slice& inner, const MX& outer) const {
    // Create return MX
    m = (*this)->get_nz_ref(ind1 ? inner-1: inner, ind1 ? outer-1.0: outer);
  }

  void MX::get_nz(MX& m, bool ind1, const MX& inner, const Slice& outer) const {
    // Create return MX
    m = (*this)->get_nz_ref(ind1 ? inner-1.0: inner, ind1 ? outer-1: outer);
  }

  void MX::set_nz(const MX& m, bool ind1, const Slice& kk) {
    // Fallback on IM
    set_nz(m, ind1, kk.all(nnz(), ind1));
  }

  void MX::set_nz(const MX& m, bool ind1, const Matrix<casadi_int>& kk) {
    casadi_assert(kk.nnz()==m.nnz() || m.nnz()==1,
      "MX::set_nz: length of non-zero indices (" + str(kk.nnz()) + ") " +
      "must match size of rhs (" + str(m.nnz()) + ").");

    // Assert dimensions of assigning matrix
    if (kk.sparsity() != m.sparsity()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        if (!m.is_dense()) return; // Nothing to set
        return set_nz(MX(kk.sparsity(), m), ind1, kk);
      } else if (kk.size() == m.size()) {
        // Project sparsity if needed
        return set_nz(project(m, kk.sparsity()), ind1, kk);
      } else if (kk.size1() == m.size2() && kk.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return set_nz(m.T(), ind1, kk);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch. lhs is " + str(kk.size())
                     + ", while rhs is " + str(m.size()));
      }
    }

    // Call recursively if points both objects point to the same node
    if (this==&m) {
      MX m_copy = m;
      return set_nz(m_copy, ind1, kk);
    }

    // Check bounds
    casadi_int sz = nnz();
    casadi_assert_in_range(kk.nonzeros(), -sz+ind1, sz+ind1);

    // Quick return if no assignments to be made
    if (kk.nnz()==0) return;

    // Handle index-1, negative indices
    if (ind1 || *std::min_element(kk->begin(), kk->end())<0) {
      Matrix<casadi_int> kk_mod = kk;
      for (auto&& i : kk_mod.nonzeros()) {
        casadi_assert(!(ind1 && i<=0),
          "Matlab is 1-based, but requested index " + str(i) + ". "
          "Note that negative slices are disabled in the Matlab interface. "
          "Possibly you may want to use 'end'.");
        if (ind1) i--;
        if (i<0) i += sz;
      }
      return set_nz(m, false, kk_mod); // Call recursively
    }

    // Create a nonzero assignment node
    *this = m->get_nzassign(*this, kk.nonzeros());
  }

  void MX::set_nz(const MX& m, bool ind1, const MX& kk) {
    *this = m->get_nzassign(*this, ind1 ? kk-1 : kk);
  }

  MX MX::binary(casadi_int op, const MX &x, const MX &y) {
    // Check, correct dimensions
    if (x.size()!=y.size() && !x.is_scalar() && !y.is_scalar()) {
      // x and y are horizontal multiples of each other?
      if (!x.is_empty() && !y.is_empty()) {
        if (x.size1() == y.size1() && x.size2() % y.size2() == 0) {
          return binary(op, x, repmat(y, 1, x.size2() / y.size2()));
        } else if (y.size1() == x.size1() && y.size2() % x.size2() == 0) {
          return binary(op, repmat(x, 1, y.size2() / x.size2()), y);
        }
      }
      // x and y are empty horizontal multiples of each other?
      if (x.size1()==0 && y.size1()==0 && x.size2()>0 && y.size2()>0) {
        if (x.size2() % y.size2() == 0) {
          return MX(0, x.size2());
        } else if (y.size2() % x.size2() == 0) {
          return MX(0, y.size2());
        }
      }
      // Dimension mismatch
      casadi_error("Dimension mismatch for " + casadi_math<double>::print(op, "x", "y") +
                   ", x is " + x.dim() + ", while y is " + y.dim());
    }
    // Call internal class
    return x->get_binary(op, y);
  }

  MX MX::unary(casadi_int op, const MX &x) {
    return x->get_unary(Operation(op));
  }

  MXNode* MX::get() const {
    return static_cast<MXNode*>(SharedObject::get());
  }

  MXNode* MX::operator->() {
    return static_cast<MXNode*>(SharedObject::operator->());
  }

  const MXNode* MX::operator->() const {
    return static_cast<const MXNode*>(SharedObject::operator->());
  }

  MX MX::inf(casadi_int nrow, casadi_int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  MX MX::inf(const std::pair<casadi_int, casadi_int> &rc) {
    return inf(rc.first, rc.second);
  }

  MX MX::inf(const Sparsity& sp) {
    return create(ConstantMX::create(sp, std::numeric_limits<double>::infinity()));
  }

  MX MX::nan(casadi_int nrow, casadi_int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  MX MX::nan(const std::pair<casadi_int, casadi_int>& rc) {
    return nan(rc.first, rc.second);
  }

  MX MX::nan(const Sparsity& sp) {
    return create(ConstantMX::create(sp, std::numeric_limits<double>::quiet_NaN()));
  }

  MX MX::eye(casadi_int n) {
    return MX(DM::eye(n));
  }

  MX MX::operator-() const {
    if ((*this)->op()==OP_NEG) {
      return (*this)->dep(0);
    } else {
      return (*this)->get_unary(OP_NEG);
    }
  }

  const Sparsity& MX::sparsity() const {
    return (*this)->sparsity();
  }

  void MX::erase(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc, bool ind1) {
    // Get sparsity of the new matrix
    Sparsity sp = sparsity();

    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sp.erase(rr, cc, ind1);

    // Create new matrix
    if (mapping.size()!=nnz()) {
      MX ret = (*this)->get_nzref(sp, mapping);
      *this = ret;
    }
  }

  void MX::erase(const std::vector<casadi_int>& rr, bool ind1) {
    // Get sparsity of the new matrix
    Sparsity sp = sparsity();

    // Erase from sparsity pattern
    std::vector<casadi_int> mapping = sp.erase(rr, ind1);

    // Create new matrix
    if (mapping.size()!=nnz()) {
      MX ret = (*this)->get_nzref(sp, mapping);
      *this = ret;
    }
  }

  void MX::enlarge(casadi_int nrow, casadi_int ncol,
                    const std::vector<casadi_int>& rr,
                    const std::vector<casadi_int>& cc, bool ind1) {
    Sparsity sp = sparsity();
    sp.enlarge(nrow, ncol, rr, cc, ind1);

    MX ret = (*this)->get_nzref(sp, range(nnz())); // FIXME?
    *this = ret;
  }

  MX MX::mtimes(const MX& x, const MX& y) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return x*y;
    } else {
      MX z = MX::zeros(Sparsity::mtimes(x.sparsity(), y.sparsity()));
      return mac(x, y, z);
    }
  }

  MX MX::einstein(const MX& A, const MX& B, const MX& C,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    return C->get_einstein(A, B, dim_c, dim_a, dim_b, c, a, b);
  }

  MX MX::einstein(const MX& A, const MX& B,
      const std::vector<casadi_int>& dim_a, const std::vector<casadi_int>& dim_b,
      const std::vector<casadi_int>& dim_c,
      const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
      const std::vector<casadi_int>& c) {
    return MX::zeros(product(dim_c), 1)->get_einstein(A, B, dim_c, dim_a, dim_b, c, a, b);
  }

  MX MX::cumsum(const MX &x, casadi_int axis) {
    if (axis==-1) axis = x.is_row();
    MX r = axis==0 ? x.T() : x;
    Sparsity sl = r(Slice(), 0).sparsity();
    MX acc = MX::sym("acc", sl);
    MX u = MX::sym("u", sl);

    Function f("f", {acc, u}, {acc+u});
    f = f.mapaccum(r.size2());
    MX ret = f(std::vector<MX>{0, r})[0];

    return axis==0 ? ret.T() : ret;
  }

  MX MX::mac(const MX& x, const MX& y, const MX& z) {
    if (x.is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + x*y;
    }

    // Check matching dimensions
    casadi_assert(x.size2()==y.size1(),
      "Matrix product with incompatible dimensions. Lhs is "
      + x.dim() + " and rhs is " + y.dim() + ".");

    // Check if we can simplify the product
    if (x.is_eye()) {
      return y + z;
    } else if (y.is_eye()) {
      return x + z;
    } else if (x.is_zero() || y.is_zero()) {
      return z;
    } else {
      return x->get_mac(y, z);
    }
  }

  MX MX::dot(const MX& x, const MX& y) {
    return x->get_dot(y);
  }

  MX MX::printme(const MX& b) const {
    return binary(OP_PRINTME, *this, b);
  }

  MX MX::attachAssert(const MX& y, const std::string &fail_message) const {
    casadi_assert(y.is_scalar(),
      "Error in attachAssert: assertion expression y must be scalar, "
      "but got " + y.dim());
    return(*this)->get_assert(y, fail_message);
  }

  MX MX::monitor(const std::string& comment) const {
    return(*this)->get_monitor(comment);
  }

  MX MX::lift(const MX& x, const MX& x_guess) {
    casadi_assert_dev(x.sparsity()==x_guess.sparsity());
    return x->_get_binary(OP_LIFT, x_guess, false, false);
  }

  DM MX::evalf(const MX& m) {
    Function f("f", std::vector<MX>{}, {m}, {{"allow_free", true}});
    return f(std::vector<DM>{})[0];
  }

  MX MX::mrdivide(const MX& b, const MX& a) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a.T(), b.T()).T();
  }

  MX MX::mldivide(const MX& a, const MX& b) {
    if (a.is_scalar() || b.is_scalar()) return b/a;
    return solve(a, b);
  }

  MX MX::dep(casadi_int ch) const {
    return (*this)->dep(ch);
  }

  casadi_int MX::n_dep() const {
    return (*this)->n_dep();
  }

  std::string MX::name() const {
    return (*this)->name();
  }

  bool MX::is_symbolic() const {
    return (*this)->op()==OP_PARAMETER;
  }

  bool MX::is_constant() const {
    return (*this)->op()==OP_CONST;
  }

  bool MX::is_call() const {
    return (*this)->op()==OP_CALL;
  }

  Function MX::which_function() const {
    return (*this)->which_function();
  }

  bool MX::is_output() const {
    return (*this)->is_output();
  }

  bool MX::has_output() const {
    return (*this)->has_output();
  }

  casadi_int MX::which_output() const {
    return (*this)->which_output();
  }

  bool MX::is_op(casadi_int op) const {
    return (*this)->op()==op;
  }

  bool MX::is_multiplication() const {
    return (*this)->op()==OP_MTIMES;
  }

  bool MX::is_norm() const {
    return dynamic_cast<const Norm*>(get())!=nullptr;
  }

  MX::operator double() const {
    return (*this)->to_double();
  }

  MX::operator DM() const {
    return (*this)->get_DM();
  }

  bool MX::is_binary() const {
    return (*this)->is_binary();
  }

  bool MX::is_unary() const {
    return (*this)->is_unary();
  }

  casadi_int MX::op() const {
    return (*this)->op();
  }

  Dict MX::info() const {
    return (*this)->info();
  }

  void MX::serialize(SerializingStream& s) const {
    return (*this)->serialize(s);
  }

  MX MX::deserialize(DeserializingStream& s) {
    return MX::create(MXNode::deserialize(s));
  }

  bool MX::is_equal(const MX& x, const MX& y, casadi_int depth) {
    return MXNode::is_equal(x.get(), y.get(), depth);
  }

  MX MX::mmin(const MX &x) {
    return x->get_mmin();
  }

  MX MX::mmax(const MX &x) {
    return x->get_mmax();
  }

  bool MX::is_commutative() const {
    if (is_unary()) return true;
    casadi_assert(is_binary() || is_unary(),
      "MX::is_commutative: must be binary or unary operation");
    return operation_checker<CommChecker>(op());
  }

  Matrix<casadi_int> MX::mapping() const {
    return (*this)->mapping();
  }

  casadi_int MX::get_temp() const {
    return (*this)->temp;
  }

  void MX::set_temp(casadi_int t) const {
    (*this)->temp = t;
  }

  casadi_int MX::n_out() const {
    return (*this)->nout();
  }

  MX MX::get_output(casadi_int oind) const {
    return (*this)->get_output(oind);
  }

  MX MX::project(const MX& x, const Sparsity& sp, bool intersect) {
    try {
      if (x.is_empty() || (sp==x.sparsity())) {
        return x;
      } else {
        casadi_assert(sp.size()==x.size(), "Cannot project " + x.dim() + " to " + sp.dim());
        if (intersect) {
          return x->get_project(sp.intersect(x.sparsity()));
        } else {
          return x->get_project(sp);
        }
      }
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("project", e.what());
    }
  }

  MX MX::densify(const MX& x, const MX& val) {
    casadi_assert_dev(val.is_scalar());
    if (x.is_dense()) {
      return x; // Already ok
    } else if (val->is_zero()) {
      return project(x, Sparsity::dense(x.size()));
    } else {
      MX ret = MX::repmat(val, x.size());
      ret(x.sparsity()) = x;
      return ret;
    }
  }

  casadi_int MX::eq_depth_ = 1;

  void MX::set_max_depth(casadi_int eq_depth) {
    eq_depth_ = eq_depth;
  }

  casadi_int MX::get_max_depth() {
    return eq_depth_;
  }

  MX MX::_sym(const std::string& name, const Sparsity& sp) {
    if (sp.nnz()==0) {
      return MX::zeros(sp);
    } else {
      return MX::create(new SymbolicMX(name, sp));
    }
  }

  bool MX::is_valid_input() const {
    return (*this)->is_valid_input();
  }

  casadi_int MX::n_primitives() const {
    return (*this)->n_primitives();
  }

  std::vector<MX> MX::primitives() const {
    std::vector<MX> ret(n_primitives());
    std::vector<MX>::iterator it=ret.begin();
    (*this)->primitives(it);
    casadi_assert_dev(it==ret.end());
    return ret;
  }

  std::vector<MX> MX::split_primitives(const MX& x) const {
    std::vector<MX> ret(n_primitives());
    std::vector<MX>::iterator it=ret.begin();
    (*this)->split_primitives(x, it);
    casadi_assert_dev(it==ret.end());
    return ret;
  }

  std::vector<SX> MX::split_primitives(const SX& x) const {
    std::vector<SX> ret(n_primitives());
    std::vector<SX>::iterator it=ret.begin();
    (*this)->split_primitives(x, it);
    casadi_assert_dev(it==ret.end());
    return ret;
  }

  std::vector<DM> MX::split_primitives(const DM& x) const {
    std::vector<DM> ret(n_primitives());
    std::vector<DM>::iterator it=ret.begin();
    (*this)->split_primitives(x, it);
    casadi_assert_dev(it==ret.end());
    return ret;
  }

  MX MX::join_primitives(const std::vector<MX>& v) const {
    casadi_assert(v.size()==n_primitives(), "Wrong number of primitives supplied");
    std::vector<MX>::const_iterator it=v.begin();
    MX ret = (*this)->join_primitives(it);
    casadi_assert_dev(it==v.end());
    return ret;
  }

  SX MX::join_primitives(const std::vector<SX>& v) const {
    casadi_assert(v.size()==n_primitives(), "Wrong number of primitives supplied");
    std::vector<SX>::const_iterator it=v.begin();
    SX ret = (*this)->join_primitives(it);
    casadi_assert_dev(it==v.end());
    return ret;
  }

  DM MX::join_primitives(const std::vector<DM>& v) const {
    casadi_assert(v.size()==n_primitives(), "Wrong number of primitives supplied");
    std::vector<DM>::const_iterator it=v.begin();
    DM ret = (*this)->join_primitives(it);
    casadi_assert_dev(it==v.end());
    return ret;
  }

  bool MX::has_duplicates() const {
    return (*this)->has_duplicates();
  }

  void MX::reset_input() const {
    (*this)->reset_input();
  }

  bool MX::is_eye() const {
    return (*this)->is_eye();
  }

  bool MX::is_zero() const {
    if (nnz()==0) {
      return true;
    } else {
      return (*this)->is_zero();
    }
  }

  bool MX::is_one() const {
    return (*this)->is_one();
  }

  bool MX::is_minus_one() const {
    return (*this)->is_value(-1);
  }

  bool MX::is_transpose() const {
    return op()==OP_TRANSPOSE;
  }

  bool MX::is_regular() const {
    if (is_constant()) {
      return static_cast<DM>(*this).is_regular();
    } else {
      casadi_error("Cannot check regularity for symbolic MX");
    }
  }

  MX MX::T() const {
    return (*this)->get_transpose();
  }

  bool MX::test_cast(const SharedObjectInternal* ptr) {
    return dynamic_cast<const MXNode*>(ptr)!=nullptr;
  }

  // Helper function
  bool has_empty(const std::vector<MX>& x, bool both=false) {
    for (auto&& i : x) {
      if (i.is_empty(both)) return true;
    }
    return false;
  }

  std::vector<MX> trim_empty(const std::vector<MX>& x, bool both=false) {
    std::vector<MX> ret;
    for (auto&& i : x) {
      if (!i.is_empty(both)) ret.push_back(i);
    }
    return ret;
  }

  MX MX::horzcat(const std::vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      std::vector<MX> ne = trim_empty(x, true);
      for (casadi_int i=0;i<ne.size();i++) {
        casadi_assert(ne[i].size1()==ne[0].size1(),
          "horzcat dimension mismatch  x[" + str(i) + "]:" + ne[i].dim() +
          " and x[0]: " + ne[0].dim() + ".");
      }
    }

    if (x.empty()) {
      return MX(1, 0);
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want horzcat(zeros(0,5),zeros(0,5)) -> zeros(0,10)
        ret = trim_empty(x, true);
        casadi_int s = 0;
        casadi_int nrow = 0;
        for (casadi_int i=0;i<ret.size();++i) {
          s+= ret[i].size2();
          casadi_assert_dev(nrow==0 || nrow==ret[i].size1());
          nrow = ret[i].size1();
        }
        return MX::zeros(nrow, s);
      } else {
        return horzcat(ret);
      }
    } else {
      return x.front()->get_horzcat(x);
    }
  }

  MX MX::diagcat(const std::vector<MX>& x) {
    if (x.empty()) {
      return MX(0, 0);
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want diagcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        casadi_int s1 = 0;
        casadi_int s2 = 0;
        for (casadi_int i=0;i<ret.size();++i) {
          s1+= ret[i].size1();
          s2+= ret[i].size2();
        }
        return MX::zeros(s1, s2);
      } else {
        return diagcat(ret);
      }
    } else {
      return x.front()->get_diagcat(x);
    }
  }

  MX MX::vertcat(const std::vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      std::vector<MX> ne = trim_empty(x, true);
      for (casadi_int i=0;i<ne.size();i++) {
        casadi_assert(ne[i].size2()==ne[0].size2(),
          "vertcat dimension mismatch  x[" + str(i) + "]:" + ne[i].dim() +
          " and x[0]: " + ne[0].dim() + ".");
      }
    }

    if (x.empty()) {
      return MX(0, 1);
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want vertcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        casadi_int s = 0;
        casadi_int ncol = 0;
        for (casadi_int i=0;i<ret.size();++i) {
          s+= ret[i].size1();
          casadi_assert_dev(ncol==0 || ret[i].size2()==ncol);
          ncol = ret[i].size2();
        }
        return MX::zeros(s, ncol);
      } else {
        return vertcat(ret);
      }
    } else if (!x.front().is_column()) {
      // Vertcat operation only supports vectors, rewrite using horzcat
      std::vector<MX> xT = x;
      for (std::vector<MX>::iterator i=xT.begin(); i!=xT.end(); ++i) *i = i->T();
      return horzcat(xT).T();
    } else {
      return x.front()->get_vertcat(x);
    }
  }

  std::vector<MX> MX::horzsplit(const MX& x, const std::vector<casadi_int>& offset) {
    // Consistency check
    casadi_assert_dev(!offset.empty());
    casadi_assert_dev(offset.front()==0);
    casadi_assert_dev(offset.back()==x.size2());
    casadi_assert_dev(is_monotone(offset));

    // Trivial return if possible
    if (offset.size()==1) {
      return std::vector<MX>(0);
    } else if (offset.size()==2) {
      return std::vector<MX>(1, x);
    } else {
      return x->get_horzsplit(offset);
    }
  }

  std::vector<MX> MX::diagsplit(const MX& x, const std::vector<casadi_int>& offset1,
                                const std::vector<casadi_int>& offset2) {
    // Consistency check
    casadi_assert_dev(!offset1.empty());
    casadi_assert_dev(offset1.front()==0);
    casadi_assert_dev(offset1.back()==x.size1());
    casadi_assert_dev(is_monotone(offset1));

    // Consistency check
    casadi_assert_dev(!offset2.empty());
    casadi_assert_dev(offset2.front()==0);
    casadi_assert_dev(offset2.back()==x.size2());
    casadi_assert_dev(is_monotone(offset2));

    return x->get_diagsplit(offset1, offset2);
  }

  std::vector<MX> MX::vertsplit(const MX& x, const std::vector<casadi_int>& offset) {
    if (x.is_column()) {
      // Consistency check
      casadi_assert_dev(!offset.empty());
      casadi_assert_dev(offset.front()==0);
      casadi_assert_dev(offset.back()==x.size1());
      casadi_assert_dev(is_monotone(offset));

      // Trivial return if possible
      if (offset.size()==1) {
        return std::vector<MX>();
      } else if (offset.size()==2) {
        return std::vector<MX>(1, x);
      } else {
        return x->get_vertsplit(offset);
      }
    } else {
      std::vector<MX> ret = horzsplit(x.T(), offset);
      for (auto&& e : ret) e = e.T();
      return ret;
    }
  }

  MX MX::blockcat(const std::vector< std::vector<MX > > &v) {
    // Quick return if no block rows
    if (v.empty()) return MX(0, 0);

    // Make sure same number of block columns
    casadi_int ncols = v.front().size();
    for (auto&& e : v) {
      casadi_assert(e.size()==ncols, "blockcat: Inconsistent number of block columns");
    }

    // Quick return if no block columns
    if (v.front().empty()) return MX(0, 0);

    // Horizontally concatenate all columns for each row, then vertically concatenate rows
    std::vector<MX> rows;
    for (auto&& e : v) {
      rows.push_back(horzcat(e));
    }
    return vertcat(rows);
  }

  MX MX::norm_2(const MX& x) {
    if (x.is_vector()) {
      return norm_fro(x);
    } else {
      return x->get_norm_2();
    }
  }

  MX MX::norm_fro(const MX& x) {
    return x->get_norm_fro();
  }

  MX MX::norm_1(const MX& x) {
    return x->get_norm_1();
  }

  MX MX::norm_inf(const MX& x) {
    return x->get_norm_inf();
  }

  MX MX::simplify(const MX& x) {
    return x;
  }

  MX MX::reshape(const MX& x, casadi_int nrow, casadi_int ncol) {
    // Quick return if trivial
    if (nrow==x.size1() && ncol==x.size2()) return x;

    // Reshape the sparsity pattern
    return reshape(x, Sparsity::reshape(x.sparsity(), nrow, ncol));
  }

  MX MX::reshape(const MX& x, const Sparsity& sp) {
    casadi_assert(sp.is_reshape(x.sparsity()), "Reshape mismatch");

    // Quick return if trivial
    if (sp==x.sparsity()) return x;

    // Call internal method
    return x->get_reshape(sp);
  }

  MX MX::sparsity_cast(const MX& x, const Sparsity& sp) {
    casadi_assert(x.nnz()==sp.nnz(),
      "Mismatching nonzero count: " + str(x.nnz()) + " versus " +
      str(sp.nnz()) + ".");

    // Quick return if trivial
    if (sp==x.sparsity()) return x;

    // Call internal method
    return x->get_sparsity_cast(sp);
  }

  MX MX::if_else(const MX &cond, const MX &x_true, const MX &x_false, bool short_circuit) {
    if (short_circuit) {
      // Get symbolic primitives
      std::vector<MX> arg = symvar(veccat(std::vector<MX>{x_true, x_false}));

      // Form functions for cases
      Function f_true("f_true", arg, {x_true});
      Function f_false("f_false", arg, {x_false});

      // Form Switch
      Function sw = Function::if_else("switch", f_true, f_false);

      // Call the Switch
      std::vector<MX> sw_arg;
      sw_arg.push_back(cond);
      sw_arg.insert(sw_arg.end(), arg.begin(), arg.end());
      return sw(sw_arg).at(0);
    } else {
      return if_else_zero(cond, x_true) + if_else_zero(!cond, x_false);
    }
  }

  MX MX::conditional(const MX& ind, const std::vector<MX>& x,
                     const MX& x_default, bool short_circuit) {
    if (short_circuit) {
      // Get symbolic primitives
      std::vector<MX> arg = x;
      arg.push_back(x_default);
      arg = symvar(veccat(arg));

      // Form functions for cases
      std::vector<Function> f(x.size());
      for (casadi_int k=0; k<x.size(); ++k) {
        std::stringstream ss;
        ss << "f_case" << k;
        f[k] = Function(ss.str(), arg, {x[k]});
      }
      Function f_default("f_default", arg, {x_default});

      // Form Switch
      Function sw = Function::conditional("switch", f, f_default);

      // Call the Switch
      std::vector<MX> sw_arg;
      sw_arg.push_back(ind);
      sw_arg.insert(sw_arg.end(), arg.begin(), arg.end());
      return sw(sw_arg).at(0);
    } else {
      MX ret = x_default;
      for (casadi_int k=0; k<x.size(); ++k) {
        ret = if_else(ind==static_cast<double>(k), x[k], ret);
      }
      return ret;
    }
  }

  MX MX::unite(const MX& A, const MX& B) {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = A.sparsity().unite(B.sparsity(), mapping);

    // Split up the mapping
    std::vector<casadi_int> nzA, nzB;

    // Copy sparsity
    for (casadi_int k=0; k<mapping.size(); ++k) {
      if (mapping[k]==1) {
        nzA.push_back(k);
      } else if (mapping[k]==2) {
        nzB.push_back(k);
      } else {
        throw CasadiException("Pattern intersection not empty");
      }
    }

    // Create mapping
    MX ret = MX::zeros(sp);
    ret = A->get_nzassign(ret, nzA);
    ret = B->get_nzassign(ret, nzB);
    return ret;
  }

  MX MX::trace(const MX& x) {
    casadi_assert(x.is_square(), "trace: must be square");
    MX res(0);
    for (casadi_int i=0; i < x.size2(); i ++) {
      res += x(i, i);
    }
    return res;
  }

  MX MX::diag(const MX& x) {
    // Nonzero mapping
    std::vector<casadi_int> mapping;

    // Get the sparsity
    Sparsity sp = x.sparsity().get_diag(mapping);

    // Create a reference to the nonzeros
    return x->get_nzref(sp, mapping);
  }

  casadi_int MX::n_nodes(const MX& x) {
    Dict opts{{"max_io", 0}, {"cse", false}, {"allow_free", true}};
    Function f("tmp_n_nodes", std::vector<MX>{}, {x}, opts);
    return f.n_nodes();
  }

  MX MX::sum2(const MX& x) {
    return mtimes(x, MX::ones(x.size2(), 1));
  }

  MX MX::sum1(const MX& x) {
    return mtimes(MX::ones(1, x.size1()), x);
  }

  MX MX::polyval(const MX& p, const MX& x) {
    casadi_assert(p.is_dense(), "polynomial coefficients vector must be a vector");
    casadi_assert(p.is_column() && p.nnz()>0, "polynomial coefficients must be a vector");
    MX ret = p.nz(0);
    for (casadi_int i=1; i<p.nnz(); ++i) {
      ret = ret*x + p.nz(i);
    }
    return ret;
  }

  std::string MX::print_operator(const MX& x, const std::vector<std::string>& args) {
    return x->disp(args);
  }

  void MX::substitute_inplace(const std::vector<MX>& v, std::vector<MX>& vdef,
                             std::vector<MX>& ex, bool reverse) {
    casadi_assert(v.size()==vdef.size(),
                          "Mismatch in the number of expression to substitute.");
    for (casadi_int k=0; k<v.size(); ++k) {
      casadi_assert(v[k].is_symbolic(),
        "Variable " + str(k) + " is not symbolic");
      casadi_assert(v[k].size() == vdef[k].size(),
        "Inconsistent shape for variable " + str(k) + ".");
    }
    casadi_assert(reverse==false, "Not implemented");

    // quick return if nothing to replace
    if (v.empty()) return;

    // implemented in MXFunction
    std::vector<MX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());
    Function temp("tmp_substitute_inplace", {v}, f_out, Dict{{"max_io", 0}, {"allow_free", true}});
    temp.get<MXFunction>()->substitute_inplace(vdef, ex);
  }

  MX MX::substitute(const MX& ex, const MX& v, const MX& vdef) {
    return substitute(std::vector<MX>{ex}, std::vector<MX>{v}, std::vector<MX>{vdef}).front();
  }

  std::vector<MX> MX::substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                                 const std::vector<MX> &vdef) {
    // Assert consistent dimensions
    casadi_assert_dev(v.size()==vdef.size());

    // Quick return if all equal
    bool all_equal = true;
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !is_equal(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Otherwise, evaluate symbolically
    Function F("tmp_substitute", v, ex, Dict{{"max_io", 0}, {"allow_free", true}});
    std::vector<MX> ret;
    F.call(vdef, ret, true);
    return ret;
  }

  MX MX::graph_substitute(const MX& x, const std::vector<MX> &v,
                          const std::vector<MX> &vdef) {
    return graph_substitute(std::vector<MX>{x}, v, vdef).at(0);
  }

  MX MX::graph_substitute(const MX& x, const std::vector<MX> &v,
                          const std::vector<MX> &vdef, bool& updated) {
    return graph_substitute(std::vector<MX>{x}, v, vdef, updated).at(0);
  }

  std::vector<MX> MX::graph_substitute(const std::vector<MX>& ex,
                                       const std::vector<MX>& v,
                                       const std::vector<MX>& vdef) {
    bool updated;
    return graph_substitute(ex, v, vdef, updated);
  }
  std::vector<MX> MX::graph_substitute(const std::vector<MX>& ex,
                                       const std::vector<MX>& v,
                                       const std::vector<MX>& vdef,
                                       bool& updated) {
    casadi_assert(v.size()==vdef.size(),
      "Mismatch in the number of expression to substitute: "
      + str(v.size()) + " <-> " + str(vdef.size()) + ".");

    updated = false;

    // Quick return if all equal
    bool all_equal = true;
    for (casadi_int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !is_equal(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Sort the expression
    Dict opts({{"max_io", 0}, {"allow_free", true}});
    Function f("tmp_graph_substitute", std::vector<MX>{}, ex, opts);
    MXFunction *ff = f.get<MXFunction>();

    // Get references to the internal data structures
    const std::vector<MXAlgEl>& algorithm = ff->algorithm_;
    std::vector<MX> swork(ff->workloc_.size()-1);

    // A boolean vector indicated whoch nodes are tainted by substitutions
    std::vector<bool> tainted(swork.size());

    // Temporary std::stringstream
    std::stringstream ss;

    // Construct lookup table for expressions,
    // giving priority to first occurances
    std::map<const MXNode*, casadi_int> expr_lookup;
    for (casadi_int i=0;i<v.size();++i) {
      auto it = expr_lookup.find(v[i].operator->());
      if (it==expr_lookup.end()) expr_lookup[v[i].operator->()] = i;
    }

    // Construct found map
    std::vector<bool> expr_found(v.size(), false);

    // Allocate output vector
    std::vector<MX> f_out(f.n_out());
    std::vector<MX> oarg, ores;

    // expr_lookup iterator
    std::map<const MXNode*, casadi_int>::const_iterator it_lookup;

    // Allocate storage for split outputs
    std::vector<std::vector<MX>> out_split(ex.size());
    for (casadi_int i = 0; i < out_split.size(); ++i) out_split[i].resize(ex[i].n_primitives());

    for (auto it=algorithm.begin(); it!=algorithm.end(); ++it) {

      if (it->op != OP_OUTPUT) {
        // Check if it->data points to a supplied expr
        it_lookup = expr_lookup.find((it->data).operator->());

        if (it_lookup!=expr_lookup.end()) {
          // Fill in that expression in-place
          MX e = vdef[it_lookup->second];

          // If node is of a MultipleOutput type
          if (e->has_output()) {
            for (casadi_int i=0;i<it->res.size();++i) {
              casadi_int k = it->res[i];
              if (k!=-1) {
                swork[k] = e.get_output(i);
                tainted[k] = true;
              }
            }
          } else {
            swork[it->res.front()] = e;
            tainted[it->res.front()] = true;
          }
          expr_found[it_lookup->second] = true;
          continue;
        } else if (it->data->has_output()) {
          bool any_tainted = false;
          // Loop over all oputputs of MultiOutput
          for (casadi_int i=0;i<it->res.size();++i) {
            // Create Output node (cached)
            casadi_int k = it->res[i];
            if (k!=-1) {
              MX out = it->data.get_output(i);
              // Check if out points to a supplied expr
              it_lookup = expr_lookup.find(out.operator->());
              if (it_lookup!=expr_lookup.end()) {
                // Fill in that expression in-place
                MX e = vdef[it_lookup->second];
                swork[k] = e;
                tainted[k] = true;
                any_tainted = true;
                expr_found[it_lookup->second] = true;
              }
            }
          }
          if (any_tainted) continue;
        }
      }

      switch (it->op) {
      case OP_INPUT:
        tainted[it->res.front()] = false;
        break;
      case OP_PARAMETER:
        swork[it->res.front()] = it->data;
        tainted[it->res.front()] = false;
        break;
      case OP_OUTPUT:
        out_split.at(it->data->ind()).at(it->data->segment()) = swork[it->arg.front()];
        break;
      default:
        {
          bool node_tainted = false;

          // Arguments of the operation
          oarg.resize(it->arg.size());
          for (casadi_int i=0; i<oarg.size(); ++i) {
            casadi_int el = it->arg[i];
            if (el>=0) node_tainted =  node_tainted || tainted[el];
            oarg[i] = el<0 ? MX(it->data->dep(i).size()) : swork.at(el);
          }

          // Perform the operation
          ores.resize(it->res.size());
          if (!node_tainted) {
            if (it->data.has_output()) {
              for (casadi_int i=0;i<it->res.size();++i) {
                ores.at(i) = it->data.get_output(i);
              }
            } else {
              ores.at(0) = it->data;
            }
          } else {
            it->data->eval_mx(oarg, ores);
          }

          // Get the result
          for (casadi_int i=0; i<ores.size(); ++i) {
            casadi_int el = it->res[i];
            if (el>=0) swork.at(el) = ores[i];
            if (el>=0) tainted[el] = node_tainted;
          }
        }
      }
    }

    // Join primitives
    for (size_t k = 0; k < out_split.size(); ++k) {
      f_out[k] = ex[k].join_primitives(out_split.at(k));
    }

    bool all_found=true;
    for (casadi_int i=0;i<v.size();++i) {
      all_found = all_found && expr_found[i];
    }

    updated = any(expr_found);

    return f_out;

  }

  void MX::extract(std::vector<MX>& ex, std::vector<MX>& v,
      std::vector<MX>& vdef, const Dict& opts) {
    try {
      // Read options
      std::string v_prefix = "v_", v_suffix = "";
      bool lift_shared = true, lift_calls = false;
      casadi_int v_ind = 0;
      for (auto&& op : opts) {
        if (op.first == "prefix") {
          v_prefix = std::string(op.second);
        } else if (op.first == "suffix") {
          v_suffix = std::string(op.second);
        } else if (op.first == "lift_shared") {
          lift_shared = op.second;
        } else if (op.first == "lift_calls") {
          lift_calls = op.second;
        } else if (op.first == "offset") {
          v_ind = op.second;
        } else {
          casadi_error("No such option: " + std::string(op.first));
        }
      }
      // Sort the expression
      Function f("tmp_extract", std::vector<MX>{}, ex, Dict{{"max_io", 0}, {"allow_free", true}});
      auto *ff = f.get<MXFunction>();
      // Get references to the internal data structures
      const std::vector<MXAlgEl>& algorithm = ff->algorithm_;
      std::vector<MX> work(ff->workloc_.size()-1);
      // Count how many times an expression has been used
      std::vector<casadi_int> usecount(work.size(), 0);
      // Remember the origin of every calculation
      std::vector<std::pair<casadi_int, casadi_int> > origin(work.size(), std::make_pair(-1, -1));
      // Which evaluations to replace
      std::vector<std::pair<casadi_int, casadi_int> > replace;
      // Evaluate the algorithm to identify which evaluations to replace
      casadi_int k=0;
      for (auto it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
        // Increase usage counters
        switch (it->op) {
        case OP_CONST:
        case OP_PARAMETER:
          break;
        default: // Unary operation, binary operation or output
          for (casadi_int c=0; c<it->arg.size(); ++c) {
            // Identify nodes used more than once
            if (lift_calls && it->op == OP_CALL) {
              // If not already marked for replacing
              if (usecount.at(it->arg[c]) >= 0) {
                replace.push_back(origin.at(it->arg[c]));
                usecount.at(it->arg[c]) = -1;  // Do not replace again
              }
            } else if (lift_shared && work[it->arg[c]].op() != OP_PARAMETER
                && work[it->arg[c]].op() != OP_CONST) {
              if (usecount.at(it->arg[c]) == 0) {
                // First time node is used
                usecount.at(it->arg[c]) = 1;
              } else if (usecount.at(it->arg[c]) == 1) {
                // Second time node is used
                replace.push_back(origin.at(it->arg[c]));
                usecount.at(it->arg[c]) = -1; // Do not replace again
              }
            }
          }
        }
        // Perform the operation
        switch (it->op) {
        case OP_OUTPUT:
          break;
        case OP_CONST:
          usecount[it->res.front()] = -1; // Never extract constants
          break;
        default:
          for (casadi_int c=0; c<it->res.size(); ++c) {
            if (it->res[c]>=0) {
              work[it->res[c]] = it->data.get_output(c);
              origin[it->res[c]] = std::make_pair(k, c);
              if (lift_calls && it->op == OP_CALL) {
                // If function call, replace right away
                replace.push_back(origin.at(it->res[c]));
                usecount.at(it->res[c]) = -1; // Do not replace again
              } else {
                usecount.at(it->res[c]) = 0; // Not (yet) extracted
              }
            }
          }
          break;
        }
      }
      // New variables and definitions
      v.clear();
      v.reserve(replace.size());
      vdef.clear();
      vdef.reserve(replace.size());
      // Quick return
      if (replace.empty()) return;
      // Sort the elements to be replaced in the order of appearence in the algorithm
      sort(replace.begin(), replace.end());
      std::vector<std::pair<casadi_int, casadi_int> >::const_iterator replace_it=replace.begin();
      // Arguments for calling the atomic operations
      std::vector<MX> oarg, ores;
      // Evaluate the algorithm
      k = 0;
      for (auto it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
        switch (it->op) {
        case OP_OUTPUT:
          casadi_assert(it->data->segment()==0, "Not implemented");
          ex[it->data->ind()] = work[it->arg.front()];
          break;
        case OP_CONST:
          work[it->res.front()] = it->data;
          break;
        default:
          {
            if (it->op == OP_PARAMETER) {
              // Free parameter
              work[it->res.front()] = it->data;
            } else {
              // Arguments of the operation
              oarg.resize(it->arg.size());
              for (casadi_int i=0; i<oarg.size(); ++i) {
                casadi_int el = it->arg[i];
                oarg[i] = el<0 ? MX(it->data->dep(i).size()) : work.at(el);
              }
              // Perform the operation
              ores.resize(it->res.size());
              it->data->eval_mx(oarg, ores);
              // Get the result
              for (casadi_int i=0; i<ores.size(); ++i) {
                casadi_int el = it->res[i];
                if (el>=0) work.at(el) = ores[i];
              }
            }
            // Possibly replace results with new variables
            for (casadi_int c=0; c<it->res.size(); ++c) {
              // Output index
              casadi_int ind = it->res[c];
              // In the list of nodes for replacing?
              bool replace_node = replace_it != replace.end()
                && replace_it->first==k && replace_it->second==c;
              // Call node (introduce variable for outputs, even if unused)
              bool output_node = lift_calls && it->op == OP_CALL;
              // Skip if no reason to replace
              if (!replace_node && !output_node) continue;
              // Create a new variable
              Sparsity v_sp = it->op == OP_PARAMETER ? it->data.sparsity() : ores.at(c).sparsity();
              v.push_back(MX::sym(v_prefix + std::to_string(v_ind++) + v_suffix, v_sp));
              // Add definition of new variable
              if (ind >= 0) {
                // Replace existing call
                casadi_assert(replace_node, "Consistency check");
                // Store the result
                vdef.push_back(work[ind]);
                // Use in calculations
                work[ind] = v.back();
                // Go to the next element to be replaced
                replace_it++;
              } else {
                // New node corresponding to an output
                casadi_assert(output_node, "Consistency check");
                // Store the result
                vdef.push_back(ores.at(c));
              }
            }
          }
        }
      }
      // Ensure all nodes have been replaced
      casadi_assert(replace_it == replace.end(), "Consistency check failed");
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("extract", e.what());
    }
  }

  void MX::shared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef,
      const std::string& v_prefix, const std::string& v_suffix) {
    // Call new, more generic function
    return extract(ex, v, vdef, Dict{{"lift_shared", true}, {"lift_calls", false},
      {"prefix", v_prefix}, {"suffix", v_suffix}});
  }

  MX MX::jacobian(const MX &f, const MX &x, const Dict& opts) {
    try {
      Dict h_opts;
      Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);
      h_opts["allow_free"] = true;
      Function h("helper_jacobian_MX", {x}, {f}, h_opts);
      return h.get<MXFunction>()->jac(opts_remainder).at(0);
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("jacobian", e.what());
    }
  }

  MX MX::hessian(const MX& f, const MX& x, const Dict& opts) {
    MX g;
    return hessian(f, x, g, opts);
  }

  MX MX::hessian(const MX& f, const MX& x, MX &g, const Dict& opts) {
    try {
      Dict all_opts = opts;
      g = gradient(f, x, opts);
      if (!opts.count("symmetric")) all_opts["symmetric"] = true;
      return jacobian(g, x, all_opts);
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("hessian", e.what());
    }
  }

  std::vector<std::vector<MX> >
  MX::forward(const std::vector<MX> &ex,
              const std::vector<MX> &arg,
              const std::vector<std::vector<MX> > &v, const Dict& opts) {
    try {
      // Read options
      bool always_inline = true;
      bool never_inline = false;

      Dict h_opts;
      Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);
      h_opts["allow_free"] = true;
      for (auto&& op : opts_remainder) {
        if (op.first=="always_inline") {
          always_inline = op.second;
        } else if (op.first=="never_inline") {
          never_inline = op.second;
        } else {
          casadi_error("No such option: " + std::string(op.first));
        }
      }
      // Call internal function on a temporary object
      Function temp("forward_temp", arg, ex, h_opts);
      std::vector<std::vector<MX> > ret;
      temp->call_forward(arg, ex, v, ret, always_inline, never_inline);
      return ret;
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("forward", e.what());
    }
  }

  std::vector<std::vector<MX> >
  MX::reverse(const std::vector<MX> &ex,
              const std::vector<MX> &arg,
              const std::vector<std::vector<MX> > &v, const Dict& opts) {
    try {
      // Read options
      bool always_inline = true;
      bool never_inline = false;


      Dict h_opts;
      Dict opts_remainder = extract_from_dict(opts, "helper_options", h_opts);
      h_opts["allow_free"] = true;

      for (auto&& op : opts_remainder) {
        if (op.first=="always_inline") {
          always_inline = op.second;
        } else if (op.first=="never_inline") {
          never_inline = op.second;
        } else {
          casadi_error("No such option: " + std::string(op.first));
        }
      }
      // Call internal function on a temporary object
      Function temp("reverse_temp", arg, ex, h_opts);
      std::vector<std::vector<MX> > ret;
      temp->call_reverse(arg, ex, v, ret, always_inline, never_inline);
      return ret;
    } catch (std::exception& e) {
      CASADI_THROW_ERROR("reverse", e.what());
    }
  }

  std::vector<bool> MX::which_depends(const MX &expr, const MX &var, casadi_int order, bool tr) {
    return _which_depends(expr, var, order, tr);
  }

  Sparsity MX::jacobian_sparsity(const MX &f, const MX &x) {
    return _jacobian_sparsity(f, x);
  }

  MX MX::det(const MX& x) {
    return x->get_det();
  }

  MX MX::inv_node(const MX& x) {
    return x->get_inv();
  }

  MX MX::inv_minor(const MX& A) {
    casadi_error("Not implemented");
  }

  MX MX::inv(const MX& x, const std::string& lsolver, const Dict& dict) {
    return solve(x, MX::eye(x.size1()), lsolver, dict);
  }

  std::vector<MX> MX::symvar(const MX& x) {
    Function f("f", std::vector<MX>{}, {x}, {{"allow_free", true}});
    return f.free_mx();
  }

  MX MX::matrix_expand(const MX& e, const std::vector<MX> &boundary, const Dict &options) {
    return matrix_expand(std::vector<MX>{e}, boundary, options).at(0);
  }

  std::vector<MX> MX::matrix_expand(const std::vector<MX>& e,
                                    const std::vector<MX> &boundary,
                                    const Dict &options) {

    // Create symbols for boundary nodes
    std::vector<MX> syms(boundary.size());

    for (casadi_int i=0;i<syms.size();++i) {
      syms[i] = MX::sym("x", boundary[i].sparsity());
    }

    // Substitute symbols for boundary nodes
    std::vector<MX> ret = graph_substitute(e, boundary, syms);

    // Obtain list of dependents
    std::vector<MX> v = symvar(veccat(ret));

    // Construct an MXFunction with it
    Function f("tmp_matrix_expand", v, ret, Dict{{"max_io", 0}, {"allow_free", true}});

    // Expand to SXFunction
    Function s = f.expand("expand_" + f.name(), options);
    std::vector<MX> r;
    s.call(graph_substitute(v, syms, boundary), r, true);
    return r;
  }

  MX MX::kron(const MX& a, const MX& b) {
    const Sparsity &a_sp = a.sparsity();
    MX filler(b.size());
    std::vector< std::vector< MX > > blocks(a.size1(), std::vector< MX >(a.size2(), filler));
    for (casadi_int i=0; i<a.size1(); ++i) {
      for (casadi_int j=0; j<a.size2(); ++j) {
        casadi_int k = a_sp.get_nz(i, j);
        if (k!=-1) {
          blocks[i][j] = a.nz(k)*b;
        }
      }
    }
    return blockcat(blocks);
  }

  MX MX::repmat(const MX& x, casadi_int n, casadi_int m) {
    if (n==0 && m==0) {
      return MX();
    } else if (n==0) {
      return MX(0, x.size2()*m);
    } else if (m==0) {
      return MX(x.size1()*n, 0);
    } else if (n==1 && m==1) {
      return x;
    } else {
      return x->get_repmat(n, m);
    }
  }

  MX MX::repsum(const MX& x, casadi_int n, casadi_int m) {
    return x->get_repsum(n, m);
  }

  MX MX::solve(const MX& a, const MX& b) {
    if (a.is_triu()) {
      // A is upper triangular
      return a->get_solve_triu(b, false);
    } else if (a.is_tril()) {
      // A is lower triangular
      return a->get_solve_tril(b, false);
    } else if (a.sparsity().is_orthonormal()) {
      // A is orthonormal -> inv(A)==A.T
      MX nz = sparsity_cast(a, Sparsity::dense(a.nnz()));
      const Sparsity& Q = a.sparsity();
      return mtimes(MX(Q, 1/nz).T(), b);
    } else {
      // Fall-back to QR factorization
      return solve(a, b, "qr");
    }
  }

  MX MX::solve(const MX& a, const MX& b, const std::string& lsolver, const Dict& dict) {
    if (a.sparsity().is_orthonormal()) return solve(a, b);
    Linsol mysolver("tmp_solve", lsolver, a.sparsity(), dict);
    return mysolver.solve(a, b, false);
  }

  MX MX::pinv(const MX& A, const std::string& lsolver, const Dict& dict) {
    if (A.size1()>=A.size2()) {
      return solve(mtimes(A.T(), A), A.T(), lsolver, dict);
    } else {
      return solve(mtimes(A, A.T()), A, lsolver, dict).T();
    }
  }

  MX MX::expm_const(const MX& A, const MX& t) {
    Dict opts;
    opts["const_A"] = true;
    Function ret = expmsol("mysolver", "slicot", A.sparsity(), opts);
    return ret(std::vector<MX>{A, t})[0];
  }

  MX MX::expm(const MX& A) {
    Function ret = expmsol("mysolver", "slicot", A.sparsity());
    return ret(std::vector<MX>{A, 1})[0];
  }

  MX MX::nullspace(const MX& A) {
    SX A_sx = SX::sym("A", A.sparsity());
    Function f("nullspace", {A_sx}, {SX::nullspace(A_sx)});
    return f(A).at(0);
  }

  bool MX::depends_on(const MX &x, const MX &arg) {
    if (x.nnz()==0) return false;

    // Construct a temporary algorithm
    Function temp("tmp_depends_on", {arg}, {x}, Dict{{"max_io", 0}, {"allow_free", true}});

    // Perform a single dependency sweep
    std::vector<bvec_t> t_in(arg.nnz(), 1), t_out(x.nnz());
    temp({get_ptr(t_in)}, {get_ptr(t_out)});

    // Loop over results
    for (casadi_int i=0; i<t_out.size(); ++i) {
      if (t_out[i]) return true;
    }

    return false;
  }

  MX MX::find(const MX& x) {
    return x->get_find();
  }

  MX MX::low(const MX& v, const MX& p, const Dict& options) {
    return p->get_low(v, options);
  }

  MX MX::bspline(const MX& x,
            const DM& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts) {
    return BSpline::create(x, knots, coeffs.nonzeros(), degree, m, opts);
  }

  MX MX::bspline(const MX& x, const MX& coeffs,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const Dict& opts) {
    return BSplineParametric::create(x, coeffs, knots, degree, m, opts);
  }

  DM MX::bspline_dual(const std::vector<double>& x,
            const std::vector< std::vector<double> >& knots,
            const std::vector<casadi_int>& degree,
            const Dict& opts) {
    return BSpline::dual(x, knots, degree, opts);
  }

  MX MX::convexify(const MX& H,
            const Dict& opts) {
    return H->get_convexify(opts);
  }


  class IncrementalSerializer {
    public:

    IncrementalSerializer() : serializer(ss, {{"debug", true}}) {
    }

    std::string pack(const MX& a) {
      // Serialization goes wrong if serialized SXNodes get destroyed
      ref.push_back(a);
      if (a.is_empty()) return "";
      // First serialize may introduce unknown dependencies (e.g. sparsity)
      // and hence definitions
      // Subsequent serialization will have references instead.
      // In order to still get a match with a later common subexpression,
      // make sure that all dependencies are already defined.
      a.serialize(serializer);
      ss.str("");
      ss.clear();
      a.serialize(serializer);
      std::string ret = ss.str();
      ss.str("");
      ss.clear();
      return ret;
    }

    private:
      std::stringstream ss;
      // List of references to keep alive
      std::vector<MX> ref;
      SerializingStream serializer;
  };


  std::vector<MX> MX::cse(const std::vector<MX>& e) {
    std::vector<MX> orig = e;
    bool updated = true;
    while (updated) {

      MX c = veccat(orig);
      Function f("f", std::vector<MX>{}, orig,
        {{"live_variables", false}, {"max_io", 0}, {"cse", false}, {"allow_free", true}});
      MXFunction *ff = f.get<MXFunction>();

      // Symbolic work, non-differentiated
      std::vector<MX> swork(ff->workloc_.size()-1);

      // Allocate storage for split outputs
      std::vector<std::vector<MX> > res_split(orig.size());
      for (casadi_int i=0; i<orig.size(); ++i) res_split[i].resize(orig[i].n_primitives());

      std::vector<MX> arg1, res1;
      std::vector<MX> res(orig.size());

      std::unordered_map<std::string, MX > cache;
      IncrementalSerializer s;

      std::unordered_map<std::string, Function> function_cache;

      // Loop over computational nodes in forward order
      casadi_int alg_counter = 0;
      for (auto it=ff->algorithm_.begin(); it!=ff->algorithm_.end(); ++it, ++alg_counter) {
        if (it->op == OP_INPUT) {
          // pass
        } else if (it->op==OP_OUTPUT) {
          // Collect the results
          res_split.at(it->data->ind()).at(it->data->segment()) = swork[it->arg.front()];
        } else if (it->op==OP_PARAMETER) {
          // Fetch parameter
          MX& target = swork[it->res.front()];
          target = it->data;
          cache[s.pack(target)] = target;
        } else {

          // Arguments of the operation
          arg1.resize(it->arg.size());
          for (casadi_int i=0; i<arg1.size(); ++i) {
            casadi_int el = it->arg[i]; // index of the argument
            arg1[i] = el<0 ? MX(it->data->dep(i).size()) : swork[el];
          }

          // Perform the operation
          res1.resize(it->res.size());
          it->data->eval_mx(arg1, res1);

          // Get the result
          for (casadi_int i=0; i<res1.size(); ++i) {
            casadi_int el = it->res[i]; // index of the output

            MX& out_i = res1[i];

            // Default assumption is that out_i is not an output node
            casadi_int output_node = -1;

            if (out_i.is_output()) {
              output_node = out_i.which_output();
              // First pack/cache the parent (MultipleOutput node e.g. Call, Horzsplit)
              out_i = out_i.dep(0);

              // If we are a call node,
              if (out_i.op()==OP_CALL) {
                std::string key = out_i.which_function().serialize();
                auto itk = function_cache.find(key);
                if (itk==function_cache.end()) {
                  function_cache[key] = out_i.which_function();
                } else {
                  out_i = Call::create_call(function_cache[key], out_i->dep_);
                }
              }
            }

            while (true) {
              // Replace out_i by a cached variant if possible
              std::string key = s.pack(out_i);

              auto itk = cache.find(key);
              if (itk==cache.end()) {
                cache[key] = out_i;
              } else {
                out_i = itk->second;
                
              }

              if (output_node==-1) {
                break; // Job is done
              } else {
                // Recreate the output node on top of the parent
                out_i = out_i.get_output(output_node);
                output_node = -1;
                // Loop once more
              }
            }

            if (el>=0) swork[el] = out_i;
          }
        }
      }

      // Join split outputs
      for (casadi_int i=0; i<res.size(); ++i) res[i] = orig[i].join_primitives(res_split[i]);

      std::vector<MX> subs_from;
      std::vector<MX> subs_to;
      for (const auto& e : function_cache) {
        e.second->merge(res, subs_from, subs_to);
      }
      orig = graph_substitute(res, subs_from, subs_to, updated);
    }

    return orig;
  }

  MX register_symbol(const MX& node, std::map<MXNode*, MX>& symbol_map,
                  std::vector<MX>& symbol_v, std::vector<MX>& parametric_v) {
    // Check if a symbol is already registered
    auto it = symbol_map.find(node.get());
    if (it==symbol_map.end()) {
      // Create a symbol and register
      MX sym = MX::sym("extracted" + str(symbol_map.size()+1), node.sparsity());
      symbol_map[node.get()] = sym;

      // Make the (symbol,parametric expression) pair available
      symbol_v.push_back(sym);
      parametric_v.push_back(node);

      // Use the new symbol
      return sym;
    } else {
      // Just use the registered symbol
      return it->second;
    }
  }

  void MX::extract_parametric(const MX &expr, const MX& par,
        MX& expr_ret, MX& symbols, MX& parametric) {

    Function f("f", {par}, {expr}, {{"live_variables", false},
      {"max_io", 0}, {"allow_free", true}});
    MXFunction *ff = f.get<MXFunction>();

    // Work vector
    std::vector< MX > w(ff->workloc_.size()-1);

    // Status of the expression:
    // 0: dependant on constants only
    // 1: dependant on parameters/constants only
    // 2: dependant on non-parameters
    std::vector< char > expr_status(ff->workloc_.size()-1, 0);

    // Split up inputs analogous to symbolic primitives
    std::vector<MX> arg_split = par.split_primitives(par);

    // Allocate storage for split outputs
    std::vector<MX> res_split;
    res_split.resize(expr.n_primitives());

    // Scratch space for node inputs/outputs
    std::vector<MX > arg1, res1;

    // Map of registered symbols
    std::map<MXNode*, MX> symbol_map;

    // Flat list of registerd symbols and parametric expressions
    std::vector<MX> symbol_v, parametric_v;

    // Loop over computational nodes in forward order
    casadi_int alg_counter = 0;
    for (auto it=ff->algorithm_.begin(); it!=ff->algorithm_.end(); ++it, ++alg_counter) {
      if (it->op == OP_INPUT) {
        w[it->res.front()] = arg_split.at(it->data->segment());
        expr_status[it->res.front()] = 1;
      } else if (it->op==OP_OUTPUT) {
        MX arg = w[it->arg.front()];
        if (expr_status[it->arg.front()]==1) {
          arg = register_symbol(arg, symbol_map, symbol_v, parametric_v);
        }
        // Collect the results
        res_split.at(it->data->segment()) = arg;
      } else if (it->op==OP_CONST) {
        // Fetch constant
        w[it->res.front()] = it->data;
        expr_status[it->res.front()] = 0;
      } else if (it->op==OP_PARAMETER) {
        // Free variables
        w[it->res.front()] = it->data;
        expr_status[it->res.front()] = 2;
      } else {
        // Arguments of the operation
        arg1.resize(it->arg.size());
        for (casadi_int i=0; i<arg1.size(); ++i) {
          casadi_int el = it->arg[i]; // index of the argument
          arg1[i] = el<0 ? MX(it->data->dep(i).size()) : w[el];
        }

        // Check worst case status of inputs
        char max_status = 0;
        for (casadi_int i=0; i<arg1.size(); ++i) {
          casadi_int el = it->arg[i]; // index of the argument
          if (el>=0) {
            max_status = std::max(max_status, expr_status[it->arg[i]]);
          }
        }
        bool any_tainted = max_status==2;

        if (any_tainted) {
          // Loop over all inputs
          for (casadi_int i=0; i<arg1.size(); ++i) {
            casadi_int el = it->arg[i]; // index of the argument

            // For each parametric input being mixed into a non-parametric expression
            if (el>=0 && expr_status[el]==1) {

              arg1[i] = register_symbol(w[el], symbol_map, symbol_v, parametric_v);
            }
          }
        }

        // Perform the operation
        res1.resize(it->res.size());
        it->data->eval_mx(arg1, res1);

        // Get the result
        for (casadi_int i=0; i<res1.size(); ++i) {
          casadi_int el = it->res[i]; // index of the output
          if (el>=0) {
            w[el] = res1[i];
            // Update expression status
            expr_status[el] = max_status;
          }
        }
      }
    }

    // Join split outputs
    expr_ret = expr.join_primitives(res_split);

    symbols = veccat(symbol_v);
    parametric = veccat(parametric_v);

  }

  MX MX::stop_diff(const MX& expr, casadi_int order) {
    std::vector<MX> s = symvar(expr);
    MX x = veccat(s);
    Dict options;
    options["never_inline"] = true;

    Dict inline_options;
    inline_options["never_inline"] = false;
    inline_options["always_inline"] = true;
    Dict der_options = Dict{{"forward_options", inline_options},
      {"reverse_options", inline_options},
      {"jacobian_options", inline_options}};
    if (order==1) {
      options["is_diff_in"] = std::vector<bool>{false};
      options["is_diff_out"] = std::vector<bool>{true};
      options = combine(options, der_options);
    } else if (order==2) {
      options["der_options"] = der_options;
      options["forward_options"] = Dict{{"is_diff_in", std::vector<bool>{false, true, true} },
        {"is_diff_out", std::vector<bool>{true}}};
      options["reverse_options"] = Dict{{"is_diff_in", std::vector<bool>{false, true, true} },
        {"is_diff_out", std::vector<bool>{true}}};
      options["jacobian_options"] = Dict{{"is_diff_in", std::vector<bool>{false, true} },
        {"is_diff_out", std::vector<bool>{false}}};
    } else {
      casadi_error("stop_diff: order must be 1 or 2, got " + str(order) + ".");
    }

    Function FS("FS", {x}, {expr}, {"x"}, {"z"}, options);
    return FS(std::vector<MX>{x})[0];
  }

  MX MX::stop_diff(const MX& expr, const MX& var, casadi_int order) {
    casadi_warning("stop_diff(expr, var, order) is not well tested.");
    std::vector<MX> xv = symvar(var);
    std::vector<MX> s = symvar(expr);
    std::vector<MX> yv = difference(s, xv);

    MX x = veccat(xv);
    MX y = veccat(yv);

    Dict options;
    options["never_inline"] = true;

    Dict inline_options;
    inline_options["never_inline"] = false;
    inline_options["always_inline"] = true;
    Dict der_options = Dict{{"forward_options", inline_options},
      {"reverse_options", inline_options},
      {"jacobian_options", inline_options}};
    if (order==1) {
      options["is_diff_in"] = std::vector<bool>{false, true};
      options["is_diff_out"] = std::vector<bool>{true};
      options = combine(options, der_options);
    } else if (order==2) {
      options["der_options"] = der_options;
      options["forward_options"] = Dict{{"is_diff_in",
        std::vector<bool>{false, true, false, true, true} },
        {"is_diff_out", std::vector<bool>{true}}};
      options["reverse_options"] = Dict{{"is_diff_in",
        std::vector<bool>{false, true, false, true}},
        {"is_diff_out", std::vector<bool>{false, true}}};
      options["jacobian_options"] = Dict{{"is_diff_in", std::vector<bool>{false, true, true}},
        {"is_diff_out", std::vector<bool>{true}}};
    } else {
      casadi_error("stop_diff: order must be 1 or 2, got " + str(order) + ".");
    }

    Function FS("FS", {x, y}, {expr}, {"x", "y"}, {"z"}, options);
    return FS(std::vector<MX>{x, y})[0];
  }

  std::vector<MX> MX::difference(const std::vector<MX>& a, const std::vector<MX>& b) {
    // Create a set of MXNodes from b
    std::set<MXNode*> bs;
    for (const auto& e : b) {
      if (!e.is_null()) bs.insert(e.get());
    }
    std::vector<MX> ret;
    for (auto&& e : a) {
      // If the element is not in the set, add it to the return vector
      if (bs.find(e.get())==bs.end()) {
        ret.push_back(e);
      }
    }
    return ret;
  }

  MX interpn_G(casadi_int i, // Dimension to interpolate along
                const MX& v, // Coefficients
                const std::vector<MX>& xis, // Normalised coordinates
                const std::vector<MX>& L, const std::vector<MX>& Lp, // Lower indices
                const std::vector<casadi_int>& strides,
                const Slice& I,
                const MX& offset=0 // Offset into coefficients vector
                ) {
    if (i==0) {
      MX ret;
      v.get_nz(ret, false, offset, I);
      return ret;
    } else {
      casadi_int j = xis.size()-i;
      MX offsetL, offsetR;
      if (strides[j]==1) {
        offsetL = offset+L[j];
        offsetR = offset+Lp[j];
      } else {
        offsetL = offset+L[j]*strides[j];
        offsetR = offsetL+strides[j];
      }
      MX vl = interpn_G(i-1, v, xis, L, Lp, strides, I, offsetL);
      MX vu = interpn_G(i-1, v, xis, L, Lp, strides, I, offsetR);

      // Perform interpolation between vl and vu
      return vl + xis[j]*(vu-vl);
    }
  }

  MX MX::interpn_linear(const std::vector<MX>& x, const MX& v, const std::vector<MX>& xq,
      const Dict& opts) {

    casadi_int n_dim = x.size();
    std::vector<std::string> lookup_mode(n_dim, "auto");
    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_mode = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    casadi_assert_dev(xq.size()==n_dim);
    casadi_assert_dev(v.is_vector());

    // Extract grid dimensions
    std::vector<casadi_int> x_dims;
    for (auto e : x) x_dims.push_back(e.numel());

    // Determine multipicity of output
    casadi_int n_out = v.numel()/product(x_dims);
    casadi_assert(n_out*product(x_dims)==v.numel(),
      "Dimension mismatch: coefficients (" + str(v.numel()) + ") should be "
      "an integer multiple of product-of-dimensions (" + str(product(x_dims)) + ").");

    // Dimension check xq
    casadi_int nq = xq[0].numel();
    for (auto e : xq) {
      casadi_assert_dev(e.is_vector() && e.numel()==nq);
    }

    // Compute stride vector
    std::vector<casadi_int> strides;
    strides.push_back(n_out);
    for (auto d : x_dims) strides.push_back(strides.back()*d);

    // Pre-compute lower index and normalized coordinate
    // (Allows for more sub-expression sharing)
    std::vector<MX> xis, Ls, Lps;
    for (casadi_int i=0;i<n_dim;++i) {
      MX L = low(x[i], xq[i], {{"lookup_mode", lookup_mode[i]}});
      MX Lp = L+1;
      MX xl, xu;
      x[i].get_nz(xl, false, L);
      x[i].get_nz(xu, false, Lp);
      xis.push_back((xq[i]-xl)/(xu-xl));
      Ls.push_back(L);
      Lps.push_back(Lp);
    }

    Slice I(0, n_out);

    return interpn_G(n_dim, v, xis, Ls, Lps, strides, I);
  }

  std::vector<MX> MX::get_input(const Function& f) {
    return f.mx_in();
  }

  std::vector<MX> MX::get_free(const Function& f) {
    return f.free_mx();
  }

  MX MX::_bilin(const MX& A, const MX& x, const MX& y) {
   return A->get_bilin(x, y);
 }

 MX MX::_rank1(const MX& A, const MX& alpha, const MX& x, const MX& y) {
   return A->get_rank1(alpha, x, y);
 }

 MX MX::_logsumexp(const MX& x) {
   return x->get_logsumexp();
 }


 void MX::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
   try {
     (*this)->eval_mx(arg, res);
   } catch (std::exception& e) {
     CASADI_THROW_ERROR_OBJ("eval_mx", e.what());
   }
 }

 void MX::ad_forward(const std::vector<std::vector<MX> >& fseed,
                 std::vector<std::vector<MX> >& fsens) const {
   try {
     (*this)->ad_forward(fseed, fsens);
   } catch (std::exception& e) {
     CASADI_THROW_ERROR_OBJ("ad_forward", e.what());
   }
 }

 void MX::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                 std::vector<std::vector<MX> >& asens) const {
   try {
     (*this)->ad_reverse(aseed, asens);
   } catch (std::exception& e) {
     CASADI_THROW_ERROR_OBJ("ad_reverse", e.what());
   }
 }

#undef CASADI_THROW_ERROR
} // namespace casadi
