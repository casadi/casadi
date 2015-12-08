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


#include "mx_node.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "multiple_output.hpp"
#include "../std_vector_tools.hpp"
#include "norm.hpp"
#include "../casadi_math.hpp"
#include "../function/mx_function.hpp"

using namespace std;
namespace casadi {

  template class GenericMatrix< MX >;

  MX::~MX() {
  }

  MX::MX() {
    assignNode(ZeroByZero::getInstance());
  }

  MX::MX(MXNode* node, bool dummy1, bool dummy2, bool dummy3, bool dummy4) {
    assignNode(node);
  }

  MX MX::create(MXNode* node) {
    return MX(node, false, false, false, false);
  }

  MX::MX(double x) {
    assignNode(ConstantMX::create(Sparsity::dense(1, 1), x));
  }

  MX::MX(const Matrix<double>& x) {
    assignNode(ConstantMX::create(x));
  }

  MX::MX(const std::vector<double>& x) {
    assignNode(ConstantMX::create(DM(x)));
  }

  MX::MX(const Sparsity& sp, const MX& val) {
    if (sp.isReshape(val.sparsity())) {
      *this = reshape(val, sp);
    } else if (val.is_scalar()) {
      // Dense matrix if val dense
      if (val.is_dense()) {
        if (val.is_constant()) {
          assignNode(ConstantMX::create(sp, val.getValue()));
        } else {
          *this = val->getGetNonzeros(sp, std::vector<int>(sp.nnz(), 0));
        }
      } else {
        // Empty matrix
        assignNode(ConstantMX::create(Sparsity(sp.size()), 0));
      }
    } else {
      casadi_assert(val.is_column() && sp.nnz()==val.size1());
      *this = densify(val)->getGetNonzeros(sp, range(sp.nnz()));
    }
  }

  MX::MX(const Sparsity& sp) {
    assignNode(ConstantMX::create(sp, 1));
  }

  MX::MX(int nrow, int ncol) {
    assignNode(ConstantMX::create(Sparsity(nrow, ncol), 0));
  }

  MX::MX(const std::pair<int, int>& rc) {
    assignNode(ConstantMX::create(Sparsity(rc), 0));
  }

  MX::MX(const Sparsity& sp, int val, bool dummy) {
    assignNode(ConstantMX::create(sp, val));
  }

  MX::MX(const Sparsity& sp, double val, bool dummy) {
    assignNode(ConstantMX::create(sp, val));
  }

  std::vector<MX> MX::createMultipleOutput(MXNode* node) {
    casadi_assert(dynamic_cast<MultipleOutput*>(node)!=0);
    MX x =  MX::create(node);
    std::vector<MX> ret(x->nout());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = MX::create(new OutputNode(x, i));
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
    return get(m, ind1, rr.getAll(size1(), ind1), cc.getAll(size2(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Slice& rr, const Matrix<int>& cc) const {
    // Fall back on (IM, IM)
    get(m, ind1, rr.getAll(size1(), ind1), cc);
  }

  void MX::get(MX& m, bool ind1, const Matrix<int>& rr, const Slice& cc) const {
    // Fall back on (IM, IM)
    get(m, ind1, rr, cc.getAll(size2(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) const {
    // Make sure dense vectors
    casadi_assert_message(rr.is_dense() && rr.is_vector(),
                          "Marix::get: First index must be a dense vector");
    casadi_assert_message(cc.is_dense() && cc.is_vector(),
                          "Marix::get: Second index must be a dense vector");

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.data(), cc.data(), mapping, ind1);

    // Create return MX
    m = (*this)->getGetNonzeros(sp, mapping);
  }

  void MX::get(MX& m, bool ind1, const Slice& rr) const {
    // Fall back on IM
    get(m, ind1, rr.getAll(numel(), ind1));
  }

  void MX::get(MX& m, bool ind1, const Matrix<int>& rr) const {
    // If the indexed matrix is dense, use nonzero indexing
    if (is_dense()) {
      return getNZ(m, ind1, rr);
    }

    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && rr.is_row()) || (is_row() && rr.is_column());

    // Get the sparsity pattern - does bounds checking
    std::vector<int> mapping;
    Sparsity sp = sparsity().sub(rr.data(), tr ? rr.sparsity().T() : rr.sparsity(),
                                 mapping, ind1);

    // Create return MX
    m = (*this)->getGetNonzeros(sp, mapping);
  }

  void MX::get(MX& m, bool ind1, const Sparsity& sp) const {
    casadi_assert_message(size()==sp.size(),
                          "get(Sparsity sp): shape mismatch. This matrix has shape "
                          << size() << ", but supplied sparsity index has shape "
                          << sp.size() << ".");
    m = project(*this, sp);
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.getAll(size1(), ind1), cc.getAll(size2(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr, const Matrix<int>& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr.getAll(size1(), ind1), cc);
  }

  void MX::set(const MX& m, bool ind1, const Matrix<int>& rr, const Slice& cc) {
    // Fall back on (IM, IM)
    set(m, ind1, rr, cc.getAll(size2(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Matrix<int>& rr, const Matrix<int>& cc) {
    // Row vector rr (e.g. in MATLAB) is transposed to column vector
    if (rr.size1()==1 && rr.size2()>1) {
      return set(m, ind1, rr.T(), cc);
    }

    // Row vector cc (e.g. in MATLAB) is transposed to column vector
    if (cc.size1()==1 && cc.size2()>1) {
      return set(m, ind1, rr, cc.T());
    }

    // Make sure rr and cc are dense vectors
    casadi_assert_message(rr.is_dense() && rr.is_column(),
                          "MX::set: First index not dense vector");
    casadi_assert_message(cc.is_dense() && cc.is_column(),
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
        casadi_error("Dimension mismatch." << "lhs is " << rr.size1() << "-by-"
                     << cc.size1() << ", while rhs is " << m.size());
      }
    }

    // Dimensions
    int sz1 = size1(), sz2 = size2();

    // Report out-of-bounds
    if (!inBounds(rr.data(), -sz1+ind1, sz1+ind1)) {
      casadi_error("set[., rr, cc] out of bounds. Your rr contains "
                   << *std::min_element(rr->begin(), rr->end()) << " up to "
                   << *std::max_element(rr->begin(), rr->end())
                   << ", which is outside the range [" << -sz1+ind1 << ","<< sz1+ind1 <<  ").");
    }
    if (!inBounds(cc.data(), -sz2+ind1, sz2+ind1)) {
      casadi_error("set [., rr, cc] out of bounds. Your cc contains "
                   << *std::min_element(cc->begin(), cc->end()) << " up to "
                   << *std::max_element(cc->begin(), cc->end())
                   << ", which is outside the range [" << -sz2+ind1 << ","<< sz2+ind1 <<  ").");
    }

    // If we are assigning with something sparse, first remove existing entries
    if (!m.is_dense()) {
      erase(rr.data(), cc.data(), ind1);
    }

    // Collect all assignments
    IM el = IM::zeros(m.sparsity());
    for (int j=0; j<el.size2(); ++j) { // Loop over columns of m
      int this_j = cc.at(j) - ind1; // Corresponding column in this
      if (this_j<0) this_j += sz2;
      for (int k=el.colind(j); k<el.colind(j+1); ++k) { // Loop over rows of m
        int i = m.row(k);
        int this_i = rr.at(i) - ind1; // Corresponding row in this
        if (this_i<0) this_i += sz1;
        el.at(k) = this_i + this_j*sz1;
      }
    }
    return set(m, false, el);
  }

  void MX::set(const MX& m, bool ind1, const Slice& rr) {
    // Fall back on IM
    set(m, ind1, rr.getAll(size1(), ind1));
  }

  void MX::set(const MX& m, bool ind1, const Matrix<int>& rr) {
    // Assert dimensions of assigning matrix
    if (rr.sparsity() != m.sparsity()) {
      if (rr.size() == m.size()) {
        // Remove submatrix to be replaced
        erase(rr.data(), ind1);

        // Find the intersection between rr's and m's sparsity patterns
        Sparsity sp = rr.sparsity() * m.sparsity();

        // Project both matrices to this sparsity
        return set(project(m, sp), ind1, project(rr, sp));
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
        casadi_error("Dimension mismatch." << "lhs is " << rr.size()
                     << ", while rhs is " << m.size());
      }
    }

    // Dimensions of this
    int sz1 = size1(), sz2 = size2(), sz = nnz(), nel = numel(), rrsz = rr.nnz();

    // Quick return if nothing to set
    if (rrsz==0) return;

    // Check bounds
    if (!inBounds(rr.data(), -nel+ind1, nel+ind1)) {
      casadi_error("set[rr] out of bounds. Your rr contains "
                   << *std::min_element(rr->begin(), rr->end()) << " up to "
                   << *std::max_element(rr->begin(), rr->end())
                   << ", which is outside the range [" << -nel+ind1 << ","<< nel+ind1 <<  ").");
    }

    // Dense mode
    if (is_dense() && m.is_dense()) {
      return setNZ(m, ind1, rr);
    }

    // Construct new sparsity pattern
    std::vector<int> new_row=sparsity().get_row(), new_col=sparsity().get_col(), nz(rr.data());
    new_row.reserve(sz+rrsz);
    new_col.reserve(sz+rrsz);
    nz.reserve(rrsz);
    for (std::vector<int>::iterator i=nz.begin(); i!=nz.end(); ++i) {
      if (ind1) (*i)--;
      if (*i<0) *i += nel;
      new_row.push_back(*i % sz1);
      new_col.push_back(*i / sz1);
    }
    Sparsity sp = Sparsity::triplet(sz1, sz2, new_row, new_col);

    // If needed, update pattern
    if (sp != sparsity()) *this = project(*this, sp);

    // Find the nonzeros corresponding to rr
    sparsity().getNZ(nz);

    // Create a nonzero assignment node
    *this = simplify(m->getSetNonzeros(*this, nz));
  }

  void MX::set(const MX& m, bool ind1, const Sparsity& sp) {
    casadi_assert_message(size()==sp.size(),
                          "set(Sparsity sp): shape mismatch. This matrix has shape "
                          << size() << ", but supplied sparsity index has shape "
                          << sp.size() << ".");
    std::vector<int> ii = sp.find();
    if (m.is_scalar()) {
      (*this)(ii) = densify(m);
    } else {
      (*this)(ii) = densify(m(ii));
    }
  }

  void MX::getNZ(MX& m, bool ind1, const Slice& kk) const {
    // Fallback on IM
    getNZ(m, ind1, kk.getAll(nnz(), ind1));
  }

  void MX::getNZ(MX& m, bool ind1, const Matrix<int>& kk) const {
    // If indexed matrix was a row/column vector, make sure that the result is too
    bool tr = (is_column() && kk.is_row()) || (is_row() && kk.is_column());

    // Quick return if no entries
    if (kk.nnz()==0) {
      m = MX::zeros(tr ? kk.sparsity().T() : kk.sparsity());
      return;
    }

    // Check bounds
    int sz = nnz();
    if (!inBounds(kk.data(), -sz+ind1, sz+ind1)) {
      casadi_error("getNZ[kk] out of bounds. Your kk contains "
                   << *std::min_element(kk->begin(), kk->end()) << " up to "
                   << *std::max_element(kk->begin(), kk->end())
                   << ", which is outside the range [" << -sz+ind1 << "," << sz+ind1 <<  ").");
    }

    // Handle index-1, negative indices
    if (ind1 || *std::min_element(kk->begin(), kk->end())<0) {
      Matrix<int> kk_mod = kk;
      for (auto&& i : kk_mod.data()) {
        casadi_assert_message(!(ind1 && i<=0), "Matlab is 1-based, but requested index " <<
                              i <<  ". Note that negative slices are" <<
                              " disabled in the Matlab interface. " <<
                              "Possibly you may want to use 'end'.");
        if (ind1) i--;
        if (i<0) i += sz;
      }
      getNZ(m, false, kk_mod); // Call recursively
      return;
    }

    // Return reference to the nonzeros
    m = (*this)->getGetNonzeros(tr ? kk.sparsity().T() : kk.sparsity(), kk.data());
  }

  void MX::setNZ(const MX& m, bool ind1, const Slice& kk) {
    // Fallback on IM
    setNZ(m, ind1, kk.getAll(nnz(), ind1));
  }

  void MX::setNZ(const MX& m, bool ind1, const Matrix<int>& kk) {
    casadi_assert_message(kk.nnz()==m.nnz() || m.nnz()==1,
                          "MX::setNZ: length of non-zero indices (" << kk.nnz() << ") " <<
                          "must match size of rhs (" << m.nnz() << ").");

    // Assert dimensions of assigning matrix
    if (kk.sparsity() != m.sparsity()) {
      if (m.is_scalar()) {
        // m scalar means "set all"
        if (!m.is_dense()) return; // Nothing to set
        return setNZ(MX(kk.sparsity(), m), ind1, kk);
      } else if (kk.size() == m.size()) {
        // Project sparsity if needed
        return setNZ(project(m, kk.sparsity()), ind1, kk);
      } else if (kk.size1() == m.size2() && kk.size2() == m.size1()
                 && std::min(m.size1(), m.size2()) == 1) {
        // m is transposed if necessary
        return setNZ(m.T(), ind1, kk);
      } else {
        // Error otherwise
        casadi_error("Dimension mismatch." << "lhs is " << kk.size()
                     << ", while rhs is " << m.size());
      }
    }

    // Call recursively if points both objects point to the same node
    if (this==&m) {
      MX m_copy = m;
      return setNZ(m_copy, ind1, kk);
    }

    // Check bounds
    int sz = nnz();
    if (!inBounds(kk.data(), -sz+ind1, sz+ind1)) {
      casadi_error("setNZ[kk] out of bounds. Your kk contains "
                   << *std::min_element(kk->begin(), kk->end()) << " up to "
                   << *std::max_element(kk->begin(), kk->end())
                   << ", which is outside the range [" << -sz+ind1 << ","<< sz+ind1 <<  ").");
    }

    // Quick return if no assignments to be made
    if (kk.nnz()==0) return;

    // Handle index-1, negative indices
    if (ind1 || *std::min_element(kk->begin(), kk->end())<0) {
      Matrix<int> kk_mod = kk;
      for (auto&& i : kk_mod.data()) {
        casadi_assert_message(!(ind1 && i<=0), "Matlab is 1-based, but requested index " <<
                              i <<  ". Note that negative slices are" <<
                              " disabled in the Matlab interface. " <<
                              "Possibly you may want to use 'end'.");
        if (ind1) i--;
        if (i<0) i += sz;
      }
      return setNZ(m, false, kk_mod); // Call recursively
    }

    // Create a nonzero assignment node
    *this = simplify(m->getSetNonzeros(*this, kk.data()));
  }

  const MX MX::at(int k) const {
    MX m;
    getNZ(m, false, k);
    return m;
  }

  /// Access a non-zero element
  NonZeros<MX, int> MX::at(int k) {
    return NonZeros<MX, int>(*this, k);
  }

  MX MX::binary(int op, const MX &x, const MX &y) {
    return x->getBinarySwitch(op, y);
  }

  MX MX::unary(int op, const MX &x) {
    return x->getUnary(Operation(op));
  }

  MXNode* MX::operator->() {
    return static_cast<MXNode*>(SharedObject::operator->());
  }

  const MXNode* MX::operator->() const {
    return static_cast<const MXNode*>(SharedObject::operator->());
  }

  MX MX::inf(int nrow, int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  MX MX::inf(const std::pair<int, int> &rc) {
    return inf(rc.first, rc.second);
  }

  MX MX::inf(const Sparsity& sp) {
    return create(ConstantMX::create(sp, numeric_limits<double>::infinity()));
  }

  MX MX::nan(int nrow, int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  MX MX::nan(const std::pair<int, int>& rc) {
    return nan(rc.first, rc.second);
  }

  MX MX::nan(const Sparsity& sp) {
    return create(ConstantMX::create(sp, numeric_limits<double>::quiet_NaN()));
  }

  MX MX::eye(int n) {
    return MX(Matrix<double>::eye(n));
  }

  MX MX::operator-() const {
    if ((*this)->op()==OP_NEG) {
      return (*this)->dep(0);
    } else {
      return (*this)->getUnary(OP_NEG);
    }
  }

  MX::MX(const MX& x) : SharedObject(x) {
  }

  const Sparsity& MX::sparsity() const {
    return (*this)->sparsity();
  }

  void MX::erase(const std::vector<int>& rr, const std::vector<int>& cc, bool ind1) {
    // Get sparsity of the new matrix
    Sparsity sp = sparsity();

    // Erase from sparsity pattern
    std::vector<int> mapping = sp.erase(rr, cc, ind1);

    // Create new matrix
    if (mapping.size()!=nnz()) {
      MX ret = (*this)->getGetNonzeros(sp, mapping);
      *this = ret;
    }
  }

  void MX::erase(const std::vector<int>& rr, bool ind1) {
    // Get sparsity of the new matrix
    Sparsity sp = sparsity();

    // Erase from sparsity pattern
    std::vector<int> mapping = sp.erase(rr, ind1);

    // Create new matrix
    if (mapping.size()!=nnz()) {
      MX ret = (*this)->getGetNonzeros(sp, mapping);
      *this = ret;
    }
  }

  void MX::enlarge(int nrow, int ncol,
                   const std::vector<int>& rr, const std::vector<int>& cc, bool ind1) {
    Sparsity sp = sparsity();
    sp.enlarge(nrow, ncol, rr, cc, ind1);

    MX ret = (*this)->getGetNonzeros(sp, range(nnz())); // FIXME?
    *this = ret;
  }

  MX MX::zz_mtimes(const MX& y) const {
    if (is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return *this*y;
    } else {
      MX z = MX::zeros(sparsity().zz_mtimes(y.sparsity()));
      return mac(*this, y, z);
    }
  }

  MX MX::zz_mac(const MX& y, const MX& z) const {
    if (is_scalar() || y.is_scalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + *this*y;
    }

    // Check matching dimensions
    casadi_assert_message(size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          << dim() << " and rhs is " << y.dim() << ".");

    // Check if we can simplify the product
    if (is_identity()) {
      return y + z;
    } else if (y.is_identity()) {
      return *this + z;
    } else if (is_zero() || y.is_zero()) {
      return z;
    } else {
      return (*this)->getMultiplication(y, z);
    }
  }

  MX MX::dot(const MX& x, const MX& y) {
    return x->getDot(y);
  }

  MX MX::zz_outer_prod(const MX& y) const {
    return mul(*this, y.T());
  }

  MX MX::zz_power(const MX& n) const {
    if (n->op()==OP_CONST) {
      return MX::binary(OP_CONSTPOW, *this, n);
    } else {
      return MX::binary(OP_POW, *this, n);
    }
  }

  MX MX::zz_constpow(const MX& b) const {
    return binary(OP_CONSTPOW, *this, b);
  }

  MX MX::zz_min(const MX& b) const {
    return binary(OP_FMIN, *this, b);
  }

  MX MX::zz_max(const MX& b) const {
    return binary(OP_FMAX, *this, b);
  }

  MX MX::zz_mod(const MX& b) const {
    return binary(OP_FMOD, *this, b);
  }

  MX MX::zz_atan2(const MX& b) const {
    return binary(OP_ATAN2, *this, b);
  }

  MX MX::printme(const MX& b) const {
    return binary(OP_PRINTME, *this, b);
  }

  MX MX::attachAssert(const MX& y, const std::string &fail_message) const {
    casadi_assert_message(y.is_scalar(),
                          "Error in attachAssert: assertion expression y must be scalar, "
                          "but got " << y.dim());
    return(*this)->getAssertion(y, fail_message);
  }

  MX MX::monitor(const std::string& comment) const {
    return(*this)->getMonitor(comment);
  }

  MX MX::zz_exp() const {
    return (*this)->getUnary(OP_EXP);
  }

  MX MX::zz_log() const {
    return (*this)->getUnary(OP_LOG);
  }

  MX MX::zz_log10() const {
    return log(*this)*(1/std::log(10.));
  }

  MX MX::zz_sqrt() const {
    return (*this)->getUnary(OP_SQRT);
  }

  MX MX::zz_sin() const {
    return (*this)->getUnary(OP_SIN);
  }

  MX MX::zz_cos() const {
    return (*this)->getUnary(OP_COS);
  }

  MX MX::zz_tan() const {
    return (*this)->getUnary(OP_TAN);
  }

  MX MX::zz_asin() const {
    return (*this)->getUnary(OP_ASIN);
  }

  MX MX::zz_acos() const {
    return (*this)->getUnary(OP_ACOS);
  }

  MX MX::zz_atan() const {
    return (*this)->getUnary(OP_ATAN);
  }

  MX MX::zz_sinh() const {
    return (*this)->getUnary(OP_SINH);
  }

  MX MX::zz_cosh() const {
    return (*this)->getUnary(OP_COSH);
  }

  MX MX::zz_tanh() const {
    return (*this)->getUnary(OP_TANH);
  }

  MX MX::zz_asinh() const {
    return (*this)->getUnary(OP_ASINH);
  }

  MX MX::zz_acosh() const {
    return (*this)->getUnary(OP_ACOSH);
  }

  MX MX::zz_atanh() const {
    return (*this)->getUnary(OP_ATANH);
  }

  MX MX::zz_floor() const {
    return (*this)->getUnary(OP_FLOOR);
  }

  MX MX::zz_ceil() const {
    return (*this)->getUnary(OP_CEIL);
  }

  MX MX::zz_abs() const {
    return (*this)->getUnary(OP_FABS);
  }

  MX MX::zz_sign() const {
    return (*this)->getUnary(OP_SIGN);
  }

  MX MX::zz_copysign(const MX& y) const {
    return MX::binary(OP_COPYSIGN, *this, y);
  }

  MX MX::zz_erfinv() const {
    return (*this)->getUnary(OP_ERFINV);
  }

  MX MX::zz_erf() const {
    return (*this)->getUnary(OP_ERF);
  }

  MX MX::zz_not() const {
    return (*this)->getUnary(OP_NOT);
  }

  MX MX::zz_lift(const MX& x_guess) const {
    casadi_assert(sparsity()==x_guess.sparsity());
    return (*this)->getBinary(OP_LIFT, x_guess, false, false);
  }

  MX MX::zz_plus(const MX& y) const {
    return MX::binary(OP_ADD, *this, y);
  }

  MX MX::zz_minus(const MX& y) const {
    return MX::binary(OP_SUB, *this, y);
  }

  MX MX::zz_times(const MX& y) const {
    return MX::binary(OP_MUL, *this, y);
  }

  MX MX::zz_rdivide(const MX& y) const {
    return MX::binary(OP_DIV, *this, y);
  }

  MX MX::zz_lt(const MX& y) const {
    return MX::binary(OP_LT, *this, y);
  }

  MX MX::zz_le(const MX& y) const {
    return MX::binary(OP_LE, *this, y);
  }

  MX MX::zz_eq(const MX& y) const {
    return MX::binary(OP_EQ, *this, y);
  }

  MX MX::zz_ne(const MX& y) const {
    return MX::binary(OP_NE, *this, y);
  }

  MX MX::zz_and(const MX& y) const {
    return MX::binary(OP_AND, *this, y);
  }

  MX MX::zz_or(const MX& y) const {
    return MX::binary(OP_OR, *this, y);
  }

  MX MX::zz_if_else_zero(const MX& y) const {
    return MX::binary(OP_IF_ELSE_ZERO, *this, y);
  }

  MX MX::zz_mrdivide(const MX& b) const {
    casadi_assert_message(is_scalar() || b.is_scalar(), "Not implemented");
    return *this/b;
  }

  MX MX::zz_mldivide(const MX& b) const {
    casadi_assert_message(is_scalar() || b.is_scalar(), "Not implemented");
    return b/ *this;
  }

  MX MX::zz_mpower(const MX& b) const {
    casadi_assert_message(is_scalar() || b.is_scalar(), "Not implemented");
    return pow(*this, b);
  }

  MX MX::getDep(int ch) const { return (*this)->dep(ch); }

  int MX::getNdeps() const { return (*this)->ndep(); }

  std::string MX::getName() const { return (*this)->getName(); }

  bool         MX::is_symbolic () const { return (*this)->op()==OP_PARAMETER; }
  bool         MX::is_constant () const { return (*this)->op()==OP_CONST; }
  bool         MX::is_call () const { return (*this)->op()==OP_CALL; }
  bool         MX::is_output () const { return (*this)->isOutputNode(); }
  int         MX::get_output () const { return (*this)->getFunctionOutput(); }
  bool         MX::is_op (int op) const { return (*this)->op()==op; }
  bool         MX::is_multiplication () const { return (*this)->op()==OP_MATMUL; }
  bool         MX::is_norm () const { return dynamic_cast<const Norm*>(get())!=0; }

  int MX::numFunctions() const { return (*this)->numFunctions(); }
  Function MX::getFunction (int i) {  return (*this)->getFunction(i); }

  double MX::getValue() const {
    return (*this)->getValue();
  }

  Matrix<double> MX::getMatrixValue() const {
    return (*this)->getMatrixValue();
  }

  bool MX::is_binary() const { return (*this)->is_binaryOp();}

  bool MX::is_unary() const { return (*this)->is_unaryOp();}

  int MX::op() const {
    return (*this)->op();
  }

  bool MX::zz_is_equal(const MX& y, int depth) const {
    return zz_is_equal(static_cast<const MXNode*>(y.get()), depth);
  }

  bool MX::zz_is_equal(const MXNode* y, int depth) const {
    if (get()==y)
      return true;
    else if (depth>0)
      return (*this)->zz_is_equal(y, depth);
    else
      return false;
  }

  bool MX::is_commutative() const {
    if (is_unary()) return true;
    casadi_assert_message(is_binary() || is_unary(),
                          "MX::is_commutative: must be binary or unary operation");
    return operation_checker<CommChecker>(op());
  }

  Matrix<int> MX::mapping() const {
    return (*this)->mapping();
  }

  int MX::getTemp() const {
    return (*this)->temp;
  }

  void MX::setTemp(int t) {
    (*this)->temp = t;
  }

  int MX::n_out() const {
    return (*this)->nout();
  }

  MX MX::getOutput(int oind) const {
    return (*this)->getOutput(oind);
  }

  MX MX::zz_project(const Sparsity& sp, bool intersect) const {
    if (is_empty() || (sp==sparsity())) {
      return *this;
    } else {
      if (intersect) {
        return (*this)->getProject(sp.intersect(sparsity()));
      } else {
        return (*this)->getProject(sp);
      }
    }
  }

  MX MX::zz_densify(const MX& val) const {
    casadi_assert(val.is_scalar());
    if (is_dense()) {
      return *this; // Already ok
    } else if (val->is_zero()) {
      return project(*this, Sparsity::dense(size()));
    } else {
      MX ret = repmat(val, size());
      ret(sparsity()) = *this;
      return ret;
    }
  }

  int MX::eq_depth_ = 1;

  void MX::setEqualityCheckingDepth(int eq_depth) {
    eq_depth_ = eq_depth;
  }

  int MX::getEqualityCheckingDepth() {
    return eq_depth_;
  }

  template<>
  MX GenericMatrix<MX>::sym(const std::string& name, const Sparsity& sp) {
    return MX::create(new SymbolicMX(name, sp));
  }

  bool MX::is_valid_input() const {
    return (*this)->is_valid_input();
  }

  int MX::n_primitives() const {
    casadi_assert_message(is_valid_input(), "Not a valid input expression");
    return (*this)->n_primitives();
  }

  std::vector<MX> MX::primitives() const {
    std::vector<MX> ret(n_primitives());
    std::vector<MX>::iterator it=ret.begin();
    (*this)->primitives(it);
    casadi_assert(it==ret.end());
    return ret;
  }

  std::vector<MX> MX::split_primitives(const MX& x) const {
    std::vector<MX> ret(n_primitives());
    std::vector<MX>::iterator it=ret.begin();
    (*this)->split_primitives(x, it);
    casadi_assert(it==ret.end());
    return ret;
  }

  MX MX::join_primitives(std::vector<MX>& v) const {
    casadi_assert_message(v.size()==n_primitives(), "Wrong number of primitives supplied");
    std::vector<MX>::const_iterator it=v.begin();
    MX ret = (*this)->join_primitives(it);
    casadi_assert(it==v.end());
    return ret;
  }

  bool MX::has_duplicates() {
    return (*this)->has_duplicates();
  }

  void MX::resetInput() {
    (*this)->resetInput();
  }

  bool MX::is_identity() const {
    return (*this)->is_identity();
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
    return (*this)->isValue(-1);
  }

  bool MX::is_transpose() const {
    return op()==OP_TRANSPOSE;
  }

  bool MX::is_regular() const {
    if (is_constant()) {
      return getMatrixValue().is_regular();
    } else {
      casadi_error("Cannot check regularity for symbolic MX");
    }
  }

  MX MX::T() const {
    return (*this)->getTranspose();
  }

  bool MX::test_cast(const SharedObjectNode* ptr) {
    return dynamic_cast<const MXNode*>(ptr)!=0;
  }

  // Helper function
  bool has_empty(const vector<MX>& x, bool both=false) {
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (i->is_empty(both)) return true;
    }
    return false;
  }

  vector<MX> trim_empty(const vector<MX>& x, bool both=false) {
    vector<MX> ret;
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (!i->is_empty(both)) ret.push_back(*i);
    }
    return ret;
  }

  MX MX::horzcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size1()==ne[0].size1(),
                      "horzcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dim() <<
                      " and x[0]: " << ne[0].dim() << ".");
      }
    }

    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want horzcat(zeros(0,5),zeros(0,5)) -> zeros(0,10)
        ret = trim_empty(x, true);
        int s = 0;
        for (int i=0;i<ret.size();++i) {
          s+= ret[i].size2();
        }
        return MX::zeros(0, s);
      } else {
        return horzcat(ret);
      }
    } else {
      return x.front()->getHorzcat(x);
    }
  }

  MX MX::diagcat(const vector<MX>& x) {
    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want diagcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        int s1 = 0;
        int s2 = 0;
        for (int i=0;i<ret.size();++i) {
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

  MX MX::vertcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size2()==ne[0].size2(),
                      "vertcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dim() <<
                      " and x[0]: " << ne[0].dim() << ".");
      }
    }

    if (x.empty()) {
      return MX();
    } else if (x.size()==1) {
      return x.front();
    } else if (has_empty(x)) {
      std::vector<MX> ret = trim_empty(x);
      if (ret.empty()) {
        // We still want vertcat(zeros(5,0),zeros(5,0)) -> zeros(10,0)
        ret = trim_empty(x, true);
        int s = 0;
        for (int i=0;i<ret.size();++i) {
          s+= ret[i].size1();
        }
        return MX::zeros(s, 0);
      } else {
        return vertcat(ret);
      }
    } else if (!x.front().is_column()) {
      // Vertcat operation only supports vectors, rewrite using horzcat
      vector<MX> xT = x;
      for (vector<MX>::iterator i=xT.begin(); i!=xT.end(); ++i) *i = i->T();
      return horzcat(xT).T();
    } else {
      return x.front()->getVertcat(x);
    }
  }

  std::vector<MX> MX::zz_horzsplit(const std::vector<int>& offset) const {
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert(offset.back()==size2());
    casadi_assert(isMonotone(offset));

    // Trivial return if possible
    if (offset.size()==1) {
      return vector<MX>(0);
    } else if (offset.size()==2) {
      return vector<MX>(1, *this);
    } else {
      return (*this)->getHorzsplit(offset);
    }
  }

  std::vector<MX> MX::zz_diagsplit(const std::vector<int>& offset1,
                                   const std::vector<int>& offset2) const {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert(offset1.back()==size1());
    casadi_assert(isMonotone(offset1));

    // Consistency check
    casadi_assert(offset2.size()>=1);
    casadi_assert(offset2.front()==0);
    casadi_assert(offset2.back()==size2());
    casadi_assert(isMonotone(offset2));

    return (*this)->get_diagsplit(offset1, offset2);
  }

  std::vector<MX> MX::zz_vertsplit(const std::vector<int>& offset) const {
    if (is_column()) {
      // Consistency check
      casadi_assert(offset.size()>=1);
      casadi_assert(offset.front()==0);
      casadi_assert(offset.back()==size1());
      casadi_assert(isMonotone(offset));

      // Trivial return if possible
      if (offset.size()==1) {
        return vector<MX>();
      } else if (offset.size()==2) {
        return vector<MX>(1, *this);
      } else {
        return (*this)->getVertsplit(offset);
      }
    } else {
      std::vector<MX> ret = horzsplit(T(), offset);
      for (std::vector<MX>::iterator it=ret.begin(); it!=ret.end(); ++it) *it = it->T();
      return ret;
    }
  }

  MX MX::blockcat(const std::vector< std::vector<MX > > &v) {
    // Quick return if no block rows
    if (v.empty()) return MX(0, 0);

    // Make sure same number of block columns
    int ncols = v.front().size();
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      casadi_assert_message(it->size()==ncols, "blockcat: Inconsistent number of blocl columns");
    }

    // Quick return if no block columns
    if (v.front().empty()) return MX(0, 0);

    // Horizontally concatenate all columns for each row, then vertically concatenate rows
    std::vector<MX> rows;
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      rows.push_back(horzcat(*it));
    }
    return vertcat(rows);
  }

  MX MX::zz_norm_2() const {
    if (is_column()) {
      return norm_F(*this);
    } else {
      return (*this)->getNorm2();
    }
  }

  MX MX::zz_norm_F() const {
    return (*this)->getNormF();
  }

  MX MX::zz_norm_1() const {
    return (*this)->getNorm1();
  }

  MX MX::zz_norm_inf() const {
    return (*this)->getNormInf();
  }

  MX MX::zz_simplify() const {
    MX ret = *this;
    if (!is_empty(true)) ret->simplifyMe(ret);
    return ret;
  }

  MX MX::zz_reshape(int nrow, int ncol) const {
    if (nrow==size1() && ncol==size2())
      return *this;
    else
      return reshape(*this, reshape(sparsity(), nrow, ncol));
  }

  MX MX::zz_reshape(const Sparsity& sp) const {
    return (*this)->getReshape(sp);
  }

  MX MX::zz_vecNZ() const {
    if (is_dense()) {
      return vec(*this);
    } else {
      return (*this)->getGetNonzeros(Sparsity::dense(nnz(), 1), range(nnz()));
    }
  }

  MX MX::zz_if_else(const MX &x_true, const MX &x_false, bool short_circuit) const {
    if (short_circuit) {
      // Get symbolic primitives
      vector<MX> arg = symvar(veccat(vector<MX>{x_true, x_false}));

      // Form functions for cases
      Function f_true("f_true", arg, {x_true});
      Function f_false("f_false", arg, {x_false});

      // Form Switch
      Function sw = Function::if_else("if_else", f_true, f_false);

      // Call the Switch
      vector<MX> sw_arg;
      sw_arg.push_back(*this);
      sw_arg.insert(sw_arg.end(), arg.begin(), arg.end());
      return sw(sw_arg).at(0);
    } else {
      return if_else_zero(*this, x_true) + if_else_zero(!*this, x_false);
    }
  }

  MX MX::zz_conditional(const std::vector<MX> &x, const MX &x_default, bool short_circuit) const {
    if (short_circuit) {
      // Get symbolic primitives
      std::vector<MX> arg = x;
      arg.push_back(x_default);
      arg = symvar(veccat(arg));

      // Form functions for cases
      vector<Function> f(x.size());
      for (int k=0; k<x.size(); ++k) {
        stringstream ss;
        ss << "f_case" << k;
        f[k] = Function(ss.str(), arg, {x[k]});
      }
      Function f_default("f_default", arg, {x_default});

      // Form Switch
      Function sw = Function::conditional("conditional", f, f_default);

      // Call the Switch
      vector<MX> sw_arg;
      sw_arg.push_back(*this);
      sw_arg.insert(sw_arg.end(), arg.begin(), arg.end());
      return sw(sw_arg).at(0);
    } else {
      MX ret = x_default;
      for (int k=0; k<x.size(); ++k) {
        ret = if_else(*this==k, x[k], ret);
      }
      return ret;
    }
  }

  MX MX::zz_unite(const MX& B) const {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = sparsity().unite(B.sparsity(), mapping);

    // Split up the mapping
    std::vector<int> nzA, nzB;

    // Copy sparsity
    for (int k=0; k<mapping.size(); ++k) {
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
    ret = (*this)->getSetNonzeros(ret, nzA);
    ret = B->getSetNonzeros(ret, nzB);
    return ret;
  }

  MX MX::zz_trace() const {
    casadi_assert_message(size2() == size1(), "trace: must be square");
    MX res(0);
    for (int i=0; i < size2(); i ++) {
      res+=(*this)(i, i);
    }
    return res;
  }

  MX MX::zz_diag() const {
    // Nonzero mapping
    std::vector<int> mapping;

    // Get the sparsity
    Sparsity sp = sparsity().get_diag(mapping);

    // Create a reference to the nonzeros
    return (*this)->getGetNonzeros(sp, mapping);
  }

  int MX::zz_countNodes() const {
    Function f("tmp", vector<MX>{}, {*this});
    return f.countNodes();
  }

  MX MX::zz_sumCols() const {
    return mul(*this, MX::ones(size2(), 1));
  }

  MX MX::zz_sumRows() const {
    return mul(MX::ones(1, size1()), *this);
  }

  MX MX::zz_polyval(const MX& x) const {
    casadi_assert_message(is_dense(), "polynomial coefficients vector must be a vector");
    casadi_assert_message(is_column() && nnz()>0, "polynomial coefficients must be a vector");
    MX ret = (*this)[0];
    for (int i=1; i<nnz(); ++i) {
      ret = ret*x + (*this)[i];
    }
    return ret;
  }

  std::string MX::zz_print_operator(const std::vector<std::string>& args) const {
    return (*this)->print(args);
  }

  void MX::zz_substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef,
                                std::vector<MX>& ex, bool reverse) {
    casadi_assert_message(v.size()==vdef.size(),
                          "Mismatch in the number of expression to substitute.");
    for (int k=0; k<v.size(); ++k) {
      casadi_assert_message(v[k].is_symbolic(), "Variable " << k << " is not symbolic");
      casadi_assert_message(v[k].size() == vdef[k].size(),
                            "Inconsistent shape for variable " << k << ".");
    }
    casadi_assert_message(reverse==false, "Not implemented");

    // quick return if nothing to replace
    if (v.empty()) return;

    // Function inputs
    std::vector<MX> f_in = v;

    // Function outputs
    std::vector<MX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    Function f("mapping", f_in, f_out);
    auto *ff = dynamic_cast<MXFunction *>(f.get());

    // Get references to the internal data structures
    std::vector<MXAlgEl>& algorithm = ff->algorithm_;
    vector<MX> work(ff->workloc_.size()-1);
    vector<MX> oarg, ores;

    for (vector<MXAlgEl>::iterator it=algorithm.begin(); it!=algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        work.at(it->res.front()) = vdef.at(it->arg.front());
        break;
      case OP_PARAMETER:
      case OP_CONST:
        work.at(it->res.front()) = it->data;
        break;
      case OP_OUTPUT:
        if (it->res.front()<vdef.size()) {
          vdef.at(it->res.front()) = work.at(it->arg.front());
        } else {
          ex.at(it->res.front()-vdef.size()) = work.at(it->arg.front());
        }
        break;
      default:
        {
          // Arguments of the operation
          oarg.resize(it->arg.size());
          for (int i=0; i<oarg.size(); ++i) {
            int el = it->arg[i];
            oarg[i] = el<0 ? MX(it->data->dep(i).size()) : work.at(el);
          }

          // Perform the operation
          ores.resize(it->res.size());
          it->data->eval_mx(oarg, ores);

          // Get the result
          for (int i=0; i<ores.size(); ++i) {
            int el = it->res[i];
            if (el>=0) work.at(el) = ores[i];
          }
        }
      }
    }
  }

  MX MX::zz_substitute(const MX& v, const MX& vdef) const {
    return substitute(vector<MX>(1, *this), vector<MX>(1, v), vector<MX>(1, vdef)).front();
  }

  std::vector<MX> MX::zz_substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                                    const std::vector<MX> &vdef) {
    // Assert consistent dimensions
    casadi_assert(v.size()==vdef.size());

    // Quick return if all equal
    bool all_equal = true;
    for (int k=0; k<v.size(); ++k) {
      if (v[k].size()!=vdef[k].size() || !is_equal(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Otherwise, evaluate symbolically
    Function F("tmp", v, ex);
    return F(vdef, true);
  }

  MX MX::zz_graph_substitute(const std::vector<MX> &v, const std::vector<MX> &vdef) const {
    return graph_substitute(std::vector<MX>(1, *this), v, vdef).at(0);
  }

  std::vector<MX> MX::zz_graph_substitute(const std::vector<MX> &ex,
                                          const std::vector<MX> &expr,
                                          const std::vector<MX> &exprs) {
    casadi_assert_message(expr.size()==exprs.size(),
                          "Mismatch in the number of expression to substitute: "
                          << expr.size() << " <-> " << exprs.size() << ".");

    // Sort the expression
    Function f("tmp", vector<MX>{}, ex);
    auto *ff = dynamic_cast<MXFunction *>(f.get());

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = ff->algorithm_;
    vector<MX> swork(ff->workloc_.size()-1);

    // A boolean vector indicated whoch nodes are tainted by substitutions
    vector<bool> tainted(swork.size());

    // Temporary stringstream
    stringstream ss;

    // Construct lookup table for expressions
    std::map<const MXNode*, int> expr_lookup;
    for (int i=0;i<expr.size();++i) {
      expr_lookup[expr[i].operator->()] = i;
    }

    // Construct found map
    std::vector<bool> expr_found(expr.size());

    // Allocate output vector
    vector<MX> f_out(f.n_out());
    vector<MX> oarg, ores;

    // expr_lookup iterator
    std::map<const MXNode*, int>::const_iterator it_lookup;

    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it!=algorithm.end(); ++it) {

      if (!(it->data).isNull()) {
        // Check if it->data points to a supplied expr
        it_lookup = expr_lookup.find((it->data).operator->());

        if (it->res.front()>=0 && it_lookup!=expr_lookup.end()) {
          // Fill in that expression in-place
          swork[it->res.front()] = exprs[it_lookup->second];
          tainted[it->res.front()] = true;
          expr_found[it_lookup->second] = true;
          continue;
        }
      }

      switch (it->op) {
      case OP_INPUT:
        tainted[it->res.front()] = false;
      case OP_PARAMETER:
        swork[it->res.front()] = it->data;
        tainted[it->res.front()] = false;
        break;
      case OP_OUTPUT:
        f_out[it->res.front()] = swork[it->arg.front()];
        break;
      default:
        {
          bool node_tainted = false;

          // Arguments of the operation
          oarg.resize(it->arg.size());
          for (int i=0; i<oarg.size(); ++i) {
            int el = it->arg[i];
            if (el>=0) node_tainted =  node_tainted || tainted[el];
            oarg[i] = el<0 ? MX(it->data->dep(i).size()) : swork.at(el);
          }

          // Perform the operation
          ores.resize(it->res.size());
          if (it->res.size()==1 && it->res[0]>=0 && !node_tainted) {
            ores.at(0) = it->data;
          } else {
            const_cast<MX&>(it->data)->eval_mx(oarg, ores);
          }

          // Get the result
          for (int i=0; i<ores.size(); ++i) {
            int el = it->res[i];
            if (el>=0) swork.at(el) = ores[i];
            if (el>=0) tainted[el] = node_tainted;
          }
        }
      }
    }

    bool all_found=true;
    for (int i=0;i<expr.size();++i) {
      all_found = all_found && expr_found[i];
    }

    //casadi_assert_message(all_found,
    //             "MXFunction::extractNodes(const std::vector<MX>& expr):"
    //             " failed to locate all input expr."
    //             << std::endl << "Here's a boolean list showing which ones where found: "
    //             << expr_found);

    return f_out;

  }

  void MX::zz_extractShared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef,
                            const std::string& v_prefix, const std::string& v_suffix) {

    // Sort the expression
    Function f("tmp", vector<MX>{}, ex);
    auto *ff = dynamic_cast<MXFunction *>(f.get());

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = ff->algorithm_;
    vector<MX> work(ff->workloc_.size()-1);

    // Count how many times an expression has been used
    vector<int> usecount(work.size(), 0);

    // Remember the origin of every calculation
    vector<pair<int, int> > origin(work.size(), make_pair(-1, -1));

    // Which evaluations to replace
    vector<pair<int, int> > replace;

    // Evaluate the algorithm to identify which evaluations to replace
    int k=0;
    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
      // Increase usage counters
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        break;
      default: // Unary operation, binary operation or output
        for (int c=0; c<it->arg.size(); ++c) {
          if (usecount[it->arg[c]]==0) {
            usecount[it->arg[c]]=1;
          } else if (usecount[it->arg[c]]==1) {
            replace.push_back(origin[it->arg[c]]);
            usecount[it->arg[c]]=-1; // Extracted, do not extract again
          }
        }
      }

      // Perform the operation
      switch (it->op) {
      case OP_OUTPUT:
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->res.front()] = -1; // Never extract since it is a primitive type
        break;
      default:
        for (int c=0; c<it->res.size(); ++c) {
          if (it->res[c]>=0) {
            work[it->res[c]] = it->data.getOutput(c);
            usecount[it->res[c]] = 0; // Not (yet) extracted
            origin[it->res[c]] = make_pair(k, c);
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
    vector<pair<int, int> >::const_iterator replace_it=replace.begin();

    // Name of intermediate variables
    stringstream v_name;

    // Arguments for calling the atomic operations
    vector<MX> oarg, ores;

    // Evaluate the algorithm
    k=0;
    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
      switch (it->op) {
      case OP_OUTPUT:     ex[it->res.front()] = work[it->arg.front()];      break;
      case OP_CONST:
      case OP_PARAMETER:  work[it->res.front()] = it->data; break;
      default:
        {
          // Arguments of the operation
          oarg.resize(it->arg.size());
          for (int i=0; i<oarg.size(); ++i) {
            int el = it->arg[i];
            oarg[i] = el<0 ? MX(it->data->dep(i).size()) : work.at(el);
          }

          // Perform the operation
          ores.resize(it->res.size());
          const_cast<MX&>(it->data)->eval_mx(oarg, ores);

          // Get the result
          for (int i=0; i<ores.size(); ++i) {
            int el = it->res[i];
            if (el>=0) work.at(el) = ores[i];
          }

          // Possibly replace results with new variables
          for (int c=0; c<it->res.size(); ++c) {
            int ind = it->res[c];
            if (ind>=0 && replace_it->first==k && replace_it->second==c) {
              // Store the result
              vdef.push_back(work[ind]);

              // Create a new variable
              v_name.str(string());
              v_name << v_prefix << v.size() << v_suffix;
              v.push_back(MX::sym(v_name.str()));

              // Use in calculations
              work[ind] = v.back();

              // Go to the next element to be replaced
              replace_it++;
            }
          }
        }
      }
    }
  }

  MX MX::jacobian(const MX &f, const MX &x, bool symmetric) {
    Function temp("helper_jacobian_MX", {x}, {f});
    return MX::jac(temp, 0, 0, false, symmetric);
  }

  MX MX::gradient(const MX& f, const MX& x) {
    Function temp("helper_gradient_MX", {x}, {f});
    return MX::grad(temp);
  }

  MX MX::tangent(const MX& f, const MX& x) {
    Function temp("helper_tangent_MX", {x}, {f});
    return MX::tang(temp);
  }

  MX MX::hessian(const MX& f, const MX& x) {
    MX g;
    return hessian(f, x, g);
  }

  MX MX::hessian(const MX& f, const MX& x, MX &g) {
    g = gradient(f, x);
    Function gfcn("gfcn", {x}, {g});
    return MX::jac(gfcn, 0, 0, false, true);
  }

  MX MX::zz_det() const {
    return (*this)->getDeterminant();
  }

  MX MX::zz_inv() const {
    return (*this)->getInverse();
  }

  std::vector<MX> MX::zz_symvar() const {
    Function f("f", vector<MX>{}, {*this});
    return f.free_mx();
  }

  MX MX::zz_matrix_expand(const MX& e, const std::vector<MX> &boundary, const Dict &options) {
    std::vector<MX> e_v(1, e);
    return matrix_expand(e_v, boundary, options).at(0);
  }

  std::vector<MX> MX::zz_matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary,
     const Dict &options) {

    // Create symbols for boundary nodes
    std::vector<MX> syms(boundary.size());

    for (int i=0;i<syms.size();++i) {
      syms[i] = MX::sym("x", boundary[i].sparsity());
    }

    // Substitute symbols for boundary nodes
    std::vector<MX> ret = graph_substitute(e, boundary, syms);

    // Obtain list of dependents
    std::vector<MX> v = symvar(veccat(ret));

    // Construct an MXFunction with it
    Function f("tmp", v, ret);

    // Expand to SXFunction
    Function s = f.expand("expand_" + f.name(), options);
    return s(graph_substitute(v, syms, boundary), true);
  }

  MX MX::zz_kron(const MX& b) const {
    const Sparsity &a_sp = sparsity();
    MX filler(b.size());
    std::vector< std::vector< MX > > blocks(size1(), std::vector< MX >(size2(), filler));
    for (int i=0; i<size1(); ++i) {
      for (int j=0; j<size2(); ++j) {
        int k = a_sp.getNZ(i, j);
        if (k!=-1) {
          blocks[i][j] = (*this)[k]*b;
        }
      }
    }
    return blockcat(blocks);
  }

  MX MX::zz_repmat(int n, int m) const {
    return (*this)->getRepmat(n, m);
  }

  MX MX::zz_repsum(int n, int m) const {
    return (*this)->getRepsum(n, m);
  }

  MX MX::zz_solve(const MX& b, const std::string& lsolver, const Dict& dict) const {
    Function mysolver = linsol("tmp", lsolver, sparsity(), b.size2(), dict);
    return mysolver.linsol_solve(*this, b, false);
  }

  MX MX::zz_pinv(const std::string& lsolver, const Dict& dict) const {
    if (size1()>=size2()) {
      return solve(mul(T(), *this), T(), lsolver, dict);
    } else {
      return solve(mul(*this, T()), *this, lsolver, dict).T();
    }
  }

  MX MX::zz_nullspace() const {
    SX n = SX::sym("A", sparsity());
    Function f("nullspace", {n}, {nullspace(n)});
    return f(*this).at(0);
  }

  bool MX::zz_dependsOn(const MX &arg) const {
    if (nnz()==0) return false;

    // Construct a temporary algorithm
    Function temp("tmp", {arg}, {*this});

    bvec_t* input_ =  get_bvec_t(temp.input().data());
    // Make a column with all variables active
    std::fill(input_, input_+temp.nnz_in(0), bvec_t(1));
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.nnz_out(0);++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  MX MX::zz_find() const {
    return (*this)->getFind();
  }

  std::vector<MX> MX::get_input(const Function& f) {
    return f.mx_in();
  }

  MX MX::jac(const Function& f, int iind, int oind, bool compact, bool symmetric) {
    return Function(f)->jac_mx(iind, oind, compact, symmetric);
  }

  MX MX::grad(const Function& f, int iind, int oind) {
    return Function(f)->grad_mx(iind, oind);
  }

  MX MX::tang(const Function& f, int iind, int oind) {
    return Function(f)->tang_mx(iind, oind);
  }

  MX MX::jac(const Function& f, const std::string & iname, int oind,
         bool compact, bool symmetric) {
    return jac(f, f.index_in(iname), oind, compact, symmetric);
  }

  MX MX::jac(const Function& f, int iind, const std::string& oname,
         bool compact, bool symmetric) {
    return jac(f, iind, f.index_out(oname), compact, symmetric);
  }

  MX MX::jac(const Function& f, const std::string& iname, const std::string& oname,
         bool compact, bool symmetric) {
    return jac(f, f.index_in(iname), f.index_out(oname), compact, symmetric);
  }

  MX MX::grad(const Function& f, const std::string& iname, int oind) {
    return grad(f, f.index_in(iname), oind);
  }

  MX MX::grad(const Function& f, int iind, const std::string& oname) {
    return grad(f, iind, f.index_out(oname));
  }

  MX MX::grad(const Function& f, const std::string& iname, const std::string& oname) {
    return grad(f, f.index_in(iname), f.index_out(oname));
  }

  MX MX::tang(const Function& f, const std::string& iname, int oind) {
    return tang(f, f.index_in(iname), oind);
  }

  MX MX::tang(const Function& f, int iind, const std::string& oname) {
    return tang(f, iind, f.index_out(oname));
  }

  MX MX::tang(const Function& f, const std::string& iname, const std::string& oname) {
    return tang(f, f.index_in(iname), f.index_out(oname));
  }

 MX MX::_bilin(const MX& A, const MX& x, const MX& y) {
   return A->getBilin(x, y);
 }

 MX MX::_rank1(const MX& A, const MX& alpha, const MX& x, const MX& y) {
   return A->getRank1(alpha, x, y);
 }

} // namespace casadi
