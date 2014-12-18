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


#include "mx.hpp"
#include "mx_node.hpp"
#include "mx_tools.hpp"
#include "../function/sx_function.hpp"
#include "call_function.hpp"
#include "symbolic_mx.hpp"
#include "constant_mx.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "norm.hpp"
#include "../casadi_math.hpp"
#include "../function/mx_function_internal.hpp"

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
    assignNode(ConstantMX::create(DMatrix(x)));
  }

  MX::MX(const Sparsity& sp, const MX& val) {
    if (sp.isReshape(val.sparsity())) {
      *this = reshape(val, sp);
    } else if (val.isScalar()) {
      // Dense matrix if val dense
      if (val.isDense()) {
        if (val.isConstant()) {
          assignNode(ConstantMX::create(sp, val.getValue()));
        } else {
          *this = val->getGetNonzeros(sp, std::vector<int>(sp.size(), 0));
        }
      } else {
        // Empty matrix
        assignNode(ConstantMX::create(Sparsity::sparse(sp.shape()), 0));
      }
    } else {
      casadi_assert(val.isVector() && sp.size()==val.size1());
      *this = dense(val)->getGetNonzeros(sp, range(sp.size()));
    }
  }

  MX::MX(const Sparsity& sp, int val) {
    assignNode(ConstantMX::create(sp, val));
  }

  MX::MX(const Sparsity& sp, double val) {
    assignNode(ConstantMX::create(sp, val));
  }

  std::vector<MX> MX::createMultipleOutput(MXNode* node) {
    casadi_assert(dynamic_cast<MultipleOutput*>(node)!=0);
    MX x =  MX::create(node);
    std::vector<MX> ret(x->getNumOutputs());
    for (int i=0; i<ret.size(); ++i) {
      ret[i] = MX::create(new OutputNode(x, i));
      if (ret[i].isEmpty(true)) {
        ret[i] = MX::sparse(0, 0);
      } else if (ret[i].size()==0) {
        ret[i] = MX::sparse(ret[i].shape());
      }
    }
    return ret;
  }

  bool MX::__nonzero__() const {
    return (*this)->__nonzero__();
  }

  const MX MX::sub(const std::vector<int>& j, int i) const {
    return sub(j, std::vector<int>(1, i));
  }

  const MX MX::sub(int rr, const std::vector<int>& cc) const {
    return sub(std::vector<int>(1, rr), cc);
  }

  const MX MX::sub(const std::vector<int>& rr, const std::vector<int>& cc) const {
    // Nonzero mapping from submatrix to full
    std::vector<int> mapping;

    // Get the sparsity pattern
    Sparsity sp = sparsity().sub(rr, cc, mapping);

    // Create return MX
    return (*this)->getGetNonzeros(sp, mapping);
  }

  const MX MX::sub(int rr, int cc) const {
    int ind = sparsity().elem(rr, cc);
    if (ind>=0) {
      return (*this)->getGetNonzeros(Sparsity::getScalar(), std::vector<int>(1, ind));
    } else {
      return (*this)->getGetNonzeros(Sparsity::getScalarSparse(), std::vector<int>(0));
    }
  }

  const MX MX::sub(const Matrix<int>& k, const std::vector<int>& ii) const {
    std::vector< int > rows = range(size1());
    std::vector< MX > temp;

    for (int i=0;i<ii.size();++i) {
      MX m(k.sparsity(), MX(0));
      for (int j=0;j<m.size();++j) {
        m[j] = sub(k.at(j), ii.at(i));
      }
      temp.push_back(m);
    }
    MX ret = horzcat(temp);
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const std::vector<int>& jj, const Matrix<int>& k) const {
    std::vector< int > cols = range(size2());
    std::vector< MX > temp;

    for (int j=0;j<jj.size();++j) {
      MX m(k.sparsity(), MX(0));
      for (int i=0;i<m.size();++i) {
        m[i] = sub(jj.at(j), k.at(i));
      }
      temp.push_back(m);
    }
    MX ret = vertcat(temp);
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const Matrix<int>& j, const Matrix<int>& i) const {
    casadi_assert_message(i.sparsity()==j.sparsity(),
                          "sub(Imatrix i, Imatrix j): sparsities must match. Got "
                          << i.dimString() << " and " << j.dimString() << ".");

    MX ret(i.sparsity(), MX(0));
    for (int k=0;k<i.size();++k) {
      ret[k] = sub(j.at(k), i.at(k));
    }
    simplify(ret);
    return ret;
  }

  const MX MX::sub(const Sparsity& sp, int dummy) const {
    casadi_assert_message(
      size2()==sp.size2() && size1()==sp.size1(),
      "sub(Sparsity sp): shape mismatch. This matrix has shape "
      << size2() << " x " << size1()
      << ", but supplied sparsity index has shape "
      << sp.size2() << " x " << sp.size1() << ".");
    std::vector<unsigned char> mappingc; // Mapping that will be filled by patternunion

    // Quick return if sparsity matches MX's sparsity
    if (sparsity()==sp) { return (*this); }

    sparsity().patternCombine(sp, false, true, mappingc);
    std::vector<int> nz(sp.size(), -1);

    int k_this = 0;     // Non-zero of this matrix
    int k_sp = 0;       // Non-zero of resulting matrix
    for (std::vector<unsigned char>::const_iterator i=mappingc.begin(); i!=mappingc.end(); ++i) {
      // In this matrix
      if (*i & 1) {
        if (*i & 4) {
          k_this++;
        } else {
          nz[k_sp++] = k_this++; // In both this matrix and in resulting matrix
        }
      } else if (*i &2) {
        k_sp++;
      }
    }

    MX ret = (*this)->getGetNonzeros(sp, nz);
    return ret;
  }

  void MX::setSub(const MX& m, int j, int i) {
    setSub(m, std::vector<int>(1, j), std::vector<int>(1, i));
  }

  void MX::setSub(const MX& m, const std::vector<int>& j, int i) {
    setSub(m, j, std::vector<int>(1, i));
  }

  void MX::setSub(const MX& m, int j, const std::vector<int>& i) {
    setSub(m, std::vector<int>(1, j), i);
  }

  void MX::setSub(const MX& m, const Slice& j, const Slice& i) {
    setSub(m, j.getAll(size1()), i.getAll(size2()));
  }

  void MX::setSub(const MX& m, const std::vector<int>& rr, const std::vector<int>& cc) {
    // Allow m to be a 1x1
    if (m.isDense() && m.isScalar()) {
      if (rr.size()>1 || cc.size()>1) {
        setSub(repmat(m, rr.size(), cc.size()), rr, cc);
        return;
      }
    }

    casadi_assert_message(rr.size()==m.size1(),
                          "Dimension mismatch." << "lhs is " << rr.size() << " x " << cc.size()
                          << ", while rhs is " << m.dimString());
    casadi_assert_message(cc.size()==m.size2(),
                          "Dimension mismatch." << "lhs is " << rr.size() << " x " << cc.size()
                          << ", while rhs is " << m.dimString());

    if (isDense() && m.isDense()) {
      // Dense mode
      int ld = size1(), ld_el = m.size1(); // leading dimensions
      std::vector<int> kk1, kk2;
      for (int i=0; i<cc.size(); ++i) {
        for (int j=0; j<rr.size(); ++j) {
          kk1.push_back(cc[i]*ld + rr[j]);
          kk2.push_back(i*ld_el+j);
        }
      }
      (*this)[kk1]=m[kk2];
    } else {
      // Sparse mode

      // Remove submatrix to be replaced
      erase(rr, cc);

      // Extend m to the same dimension as this
      MX el_ext = m;
      el_ext.enlarge(size1(), size2(), rr, cc);

      // Unite the sparsity patterns
      *this = unite(*this, el_ext);
    }
  }

  void MX::setSub(const MX& m, const std::vector<int>& jj, const Matrix<int>& i) {
    // If m is scalar
    if (m.isScalar() && (jj.size() > 1 || i.size() > 1)) {
      setSub(repmat(MX(i.sparsity(), m), 1, jj.size()), jj, i);
      return;
    }

    if (!inBounds(jj, size1())) {
      casadi_error(
        "setSub[., i, jj] out of bounds. Your jj contains "
        << *std::min_element(jj.begin(), jj.end()) << " up to "
        << *std::max_element(jj.begin(), jj.end())
        << ", which is outside of the matrix shape " << dimString() << ".");
    }

    //Sparsity result_sparsity = repmat(i, 1, jj.size()).sparsity();
    Sparsity result_sparsity = vertcat(std::vector< Matrix<int> >(jj.size(), i)).sparsity();

    casadi_assert_message(
      result_sparsity == m.sparsity(),
      "setSub(., Imatrix" << i.dimString() << ", Ivector(length=" << jj.size()
      << "), Matrix<T>)::Dimension mismatch. The sparsity of repmat(Imatrix,1, "
      << jj.size() << ") = " << result_sparsity.dimString()
      << " must match the sparsity of MX = "  << m.dimString() << ".");

    std::vector<int> slice_i = range(i.size2());

    for (int k=0; k<jj.size(); ++k) {
      MX el_k = m(range(k*i.size1(), (k+1)*i.size1()), slice_i);
      for (int j=0;j<i.size();++j) {
        (*this)(jj[k], i.at(j))=el_k[j];
      }
    }

  }

  void MX::setSub(const MX& m, const Matrix<int>& j, const std::vector<int>& ii) {
    // If m is scalar
    if (m.isScalar() && (ii.size() > 1 || j.size() > 1)) {
      setSub(repmat(MX(j.sparsity(), m), ii.size(), 1), j, ii);
      return;
    }

    if (!inBounds(ii, size2())) {
      casadi_error("setSub[., ii, j] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    //Sparsity result_sparsity = repmat(j, ii.size(), 1).sparsity();
    Sparsity result_sparsity = horzcat(std::vector< Matrix<int> >(ii.size(), j)).sparsity();

    casadi_assert_message(
      result_sparsity == m.sparsity(),
      "setSub(Ivector(length=" << ii.size() << "), Imatrix" << j.dimString()
      << ", MX)::Dimension mismatch. The sparsity of repmat(Imatrix, "
      << ii.size() << ",1) = " << result_sparsity.dimString()
      << " must match the sparsity of Matrix<T> = " << m.dimString() << ".");

    std::vector<int> slice_j = range(j.size1());

    for (int k=0; k<ii.size(); ++k) {
      MX el_k = m(slice_j, range(k*j.size2(), (k+1)*j.size2()));
      for (int i=0;i<j.size();++i) {
        (*this)(j.at(i), ii[k])=el_k[i];
      }
    }

  }


  void MX::setSub(const MX& m, const Matrix<int>& j, const Matrix<int>& i) {
    casadi_assert_message(
      i.sparsity()==j.sparsity(),
      "setSub(Imatrix m, Imatrix i, Imatrix j): sparsities must match. Got "
      << i.dimString() << " for i and " << j.dimString() << " for j.");

    // If m is scalar
    if (m.isScalar() && i.numel() > 1) {
      setSub(MX(i.sparsity(), m), j, i);
      return;
    }

    casadi_assert_message(
      m.sparsity()==i.sparsity(),
      "setSub(MX m, Imatrix i, Imatrix j): sparsities must match. Got "
      << m.dimString() << " for m and " << j.dimString() << " for i and j.");

    for (int k=0; k<i.size(); ++k) {
      (*this)(j.at(k), i.at(k)) = m[k];
    }
  }

  void MX::setSub(const MX& m, const Sparsity& sp, int dummy) {
    casadi_assert_message(
      size2()==sp.size2() && size1()==sp.size1(),
      "setSub(., Sparsity sp): shape mismatch. This matrix has shape "
      << size2() << " x " << size1()
      << ", but supplied sparsity index has shape "
      << sp.size2() << " x " << sp.size1() << ".");

    // If m is scalar
    if (m.isScalar()) {
      setSub(MX(sp, m), sp, dummy);
      return;
    }

    MX mm = m.sub(sp);

    std::vector<unsigned char> mappingc; // Mapping that will be filled by patternunion

    sparsity().patternCombine(sp, false, true, mappingc);
    std::vector<int> nz(sp.size(), -1);

    int k_this = 0;     // Non-zero of this matrix
    int k_sp = 0;       // Non-zero of resulting matrix
    for (std::vector<unsigned char>::const_iterator i=mappingc.begin(); i!=mappingc.end(); ++i) {
      // In this matrix
      if (*i & 1) {
        if (*i & 4) {
          k_this++;
        } else {
          nz[k_sp++] = k_this++; // In both this matrix and in resulting matrix
        }
      } else if (*i &2) {
        k_sp++;
      }
    }

    *this =  mm->getSetNonzeros((*this), nz);

  }

  MX MX::getNZ(int k) const {
    if (k<0) k+=size();
    casadi_assert_message(k<size(),
                          "MX::getNZ: requested at(" <<  k << "), but that is out of bounds:  "
                          << dimString() << ".");
    return getNZ(std::vector<int>(1, k));
  }

  MX MX::getNZ(const std::vector<int>& k) const {
    Sparsity sp = Sparsity::dense(k.size());

    for (int i=0;i<k.size();i++) {
      casadi_assert_message(k[i] < size(), "Mapping::assign: index vector reaches " << k[i]
                            << ", while dependent is only of size " << size());
    }

    MX ret = (*this)->getGetNonzeros(sp, k);
    return ret;
  }

  MX MX::getNZ(const Matrix<int>& k) const {
    MX ret = (*this)->getGetNonzeros(k.sparsity(), k.data());
    return ret;
  }

  void MX::setNZ(int k, const MX& el) {
    if (k<0) k+=size();
    casadi_assert_message(k<size(),
                          "MX::setNZ: requested at(" <<  k << "), but that is out of bounds:  "
                          << dimString() << ".");
    setNZ(std::vector<int>(1, k), el);
  }

  void MX::setNZ(const std::vector<int>& k, const MX& el) {
    casadi_assert_message(k.size()==el.size() || el.size()==1,
                          "MX::setNZ: length of non-zero indices (" << k.size() << ") " <<
                          "must match size of rhs (" << el.size() << ").");

    // Call recursively if points both objects point to the same node
    if (this==&el) {
      MX el2 = el;
      setNZ(k, el2);
      return;
    }

    // Assert correctness
    for (int i=0; i<k.size(); ++i) {
      casadi_assert_message(k[i] < size(),
                            "Mapping::assign: index vector reaches " << k[i]
                            << ", while dependent is only of size " << size());
    }

    // Quick return if no assignments to be made
    if (k.empty()) return;

    // Temporary
    MX x;

    // Project scalars
    if (k.size()!=el.size() && el.isScalar() && el.isDense()) {
      MX new_el = el->getGetNonzeros(Sparsity::dense(1, k.size()), std::vector<int>(k.size(), 0));
      x = new_el->getSetNonzeros(*this, k);
    } else {
      // Create a nonzero assignment node
      x = el->getSetNonzeros(*this, k);
    }
    *this = x;
    simplify(*this);
  }

  void MX::setNZ(const Matrix<int>& kk, const MX& m) {
    if (m.size()==1 && m.numel()==1) {
      setNZ(kk.data(), m);
      return;
    }
    casadi_assert_message(kk.sparsity()==m.sparsity(),
                          "Matrix<T>::setNZ: sparsity of IMatrix index " << kk.dimString()
                          << " " << std::endl << "must match sparsity of rhs " << m.dimString()
                          << ".");
    setNZ(kk.data(), m);
  }

  const MX MX::at(int k) const {
    return getNZ(k);
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

  MX MX::repmat(const MX& x, const Sparsity& sp) {
    casadi_assert_message(x.isScalar(), "repmat(MX x, Sparsity sp) only defined for scalar x");
    return MX(sp, x);
  }

  MX MX::repmat(const MX& x, const std::pair<int, int> &rc) {
    return repmat(x, rc.first, rc.second);
  }

  MX MX::repmat(const MX& x, int nrow, int ncol) {
    if (x.isScalar()) {
      return MX(Sparsity::dense(nrow, ncol), x);
    } else {
      casadi_assert_message(0, "not implemented");
      return MX();
    }
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
    Matrix<double> I(Sparsity::diag(n), 1);
    return MX(I);
  }


  MX MX::operator-() const {

    if ((*this)->getOp()==OP_NEG) {
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

  Sparsity& MX::sparsityRef() {
    // Since we can potentially change the behavior of the MX node,
    // we must make a deep copy if there are other references
    makeUnique();

    // Return the reference, again, deep copy if multiple references
    (*this)->sparsity_.makeUnique();
    return (*this)->sparsity_;
  }

  void MX::erase(const std::vector<int>& rr, const std::vector<int>& cc) {
    // Get sparsity of the new matrix
    Sparsity sp = sparsity();

    // Erase from sparsity pattern
    std::vector<int> mapping = sp.erase(rr, cc);

    // Create new matrix
    if (mapping.size()!=size()) {
      MX ret = (*this)->getGetNonzeros(sp, mapping);
      *this = ret;
    }
  }

  void MX::enlarge(int nrow, int ncol, const std::vector<int>& rr, const std::vector<int>& cc) {
    Sparsity sp = sparsity();
    sp.enlarge(nrow, ncol, rr, cc);

    MX ret = (*this)->getGetNonzeros(sp, range(size()));
    *this = ret;
  }

  MX MX::zz_mtimes(const MX& y) const {
    MX z = MX::zeros(sparsity().patternProduct(y.sparsity()));
    return zz_mtimes(y, z);
  }

  MX MX::zz_mtimes(const MX& y, const MX& z) const {
    if (isScalar() || y.isScalar()) {
      // Use element-wise multiplication if at least one factor scalar
      return z + *this*y;
    }

    // Check matching dimensions
    casadi_assert_message(size2()==y.size1(),
                          "Matrix product with incompatible dimensions. Lhs is "
                          << dimString() << " and rhs is " << y.dimString() << ".");

    // Check if we can simplify the product
    if (isIdentity()) {
      return y + z;
    } else if (y.isIdentity()) {
      return *this + z;
    } else if (isZero() || y.isZero()) {
      return z;
    } else {
      return (*this)->getMultiplication(y, z);
    }
  }

  MX MX::inner_prod(const MX& y) const {
    return (*this)->getInnerProd(y);
  }

  MX MX::outer_prod(const MX& y) const {
    return mul(*this, y.T());
  }

  MX MX::zz_power(const MX& n) const {
    if (n->getOp()==OP_CONST) {
      return MX::binary(OP_CONSTPOW, *this, n);
    } else {
      return MX::binary(OP_POW, *this, n);
    }
  }

  MX MX::constpow(const MX& b) const {
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
    casadi_assert_message(y.isScalar(),
                          "Error in attachAssert: assertion expression y must be scalar, "
                          "but got " << y.dimString());
    return(*this)->getAssertion(y, fail_message);
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

  MX MX::__copysign__(const MX& y) const {
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

  void MX::lift(const MX& x_guess) {
    casadi_assert(sparsity()==x_guess.sparsity());
    *this = (*this)->getBinary(OP_LIFT, x_guess, false, false);
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

  MX MX::if_else_zero(const MX& y) const {
    return MX::binary(OP_IF_ELSE_ZERO, *this, y);
  }

  MX MX::__constpow__(const MX& b) const { return (*this).constpow(b);}
  MX MX::__mrdivide__(const MX& b) const
  { if (b.isScalar()) return *this/b; throw CasadiException("mrdivide: Not implemented");}
  MX MX::zz_mpower(const MX& b) const {
    return pow(*this, b); throw CasadiException("mpower: Not implemented");
  }

  void MX::append(const MX& y) {
    *this = vertcat(*this, y);
  }

  void MX::appendColumns(const MX& y) {
    *this = horzcat(*this, y);
  }

  long MX::max_num_calls_in_print_ = 10000;

  void MX::setMaxNumCallsInPrint(long num) {
    max_num_calls_in_print_ = num;
  }

  long MX::getMaxNumCallsInPrint() {
    return max_num_calls_in_print_;
  }

  MX MX::getDep(int ch) const { return (*this)->dep(ch); }

  int MX::getNdeps() const { return (*this)->ndep(); }

  std::string MX::getName() const { return (*this)->getName(); }

  bool         MX::isSymbolic () const { return (*this)->getOp()==OP_PARAMETER; }
  bool         MX::isConstant () const { return (*this)->getOp()==OP_CONST; }
  bool         MX::isEvaluation () const { return (*this)->getOp()==OP_CALL; }
  bool         MX::isEvaluationOutput () const { return (*this)->isOutputNode(); }

  int         MX::getEvaluationOutput () const { return (*this)->getFunctionOutput(); }


  bool         MX::isOperation (int op) const { return (*this)->getOp()==op; }
  bool         MX::isMultiplication () const { return (*this)->getOp()==OP_MATMUL; }

  bool         MX::isNorm () const { return dynamic_cast<const Norm*>(get())!=0; }

  Function MX::getFunction () {  return (*this)->getFunction(); }

  double MX::getValue() const {
    return (*this)->getValue();
  }

  Matrix<double> MX::getMatrixValue() const {
    return (*this)->getMatrixValue();
  }

  bool MX::isBinary() const { return (*this)->isBinaryOp();}

  bool MX::isUnary() const { return (*this)->isUnaryOp();}

  int MX::getOp() const {
    return (*this)->getOp();
  }

  bool MX::isEqual(const MX& y, int depth) const {
    return isEqual(static_cast<const MXNode*>(y.get()), depth);
  }

  bool MX::isEqual(const MXNode* y, int depth) const {
    if (get()==y)
      return true;
    else if (depth>0)
      return (*this)->isEqual(y, depth);
    else
      return false;
  }

  bool MX::isCommutative() const {
    if (isUnary()) return true;
    casadi_assert_message(isBinary() || isUnary(),
                          "MX::isCommutative: must be binary or unary operation");
    return operation_checker<CommChecker>(getOp());
  }

  long MX::__hash__() const {
    // TODO(Joel):
    // This is bad coding, the pointer might on certain architectures be larger
    // than the long, giving an error of type "error: cast from 'const SharedObjectInternal*'
    // to 'int' loses precision". The solution is to convert to intptr_t, but this is not
    // yet C++ standard
    return long(get());
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

  int MX::getNumOutputs() const {
    return (*this)->getNumOutputs();
  }

  MX MX::getOutput(int oind) const {
    return (*this)->getOutput(oind);
  }

  MX MX::setSparse(const Sparsity& sp, bool intersect) const {
    if (isEmpty() || (sp==sparsity())) {
      return *this;
    } else {
      if (intersect) {
        return (*this)->getSetSparse(sp.patternIntersection(sparsity()));
      } else {
        return (*this)->getSetSparse(sp);
      }
    }
  }

  void MX::densify(const MX& val) {
    casadi_assert(val.isScalar());
    if (isDense()) {
      return; // Already ok
    } else if (val->isZero()) {
      *this = setSparse(Sparsity::dense(shape()));
    } else {
      MX new_this = repmat(val, shape());
      new_this(sparsity()) = *this;
      *this = new_this;
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

  bool MX::isSymbolicSparse() const {
    if (getOp()==OP_HORZCAT) {
      // Check if the expression is a horzcat where all components are symbolic primitives
      for (int d=0; d<getNdeps(); ++d) {
        if (!getDep(d).isSymbolic()) {
          return false;
        }
      }
      return true;
    } else {
      return isSymbolic();
    }
  }

  bool MX::isIdentity() const {
    return (*this)->isIdentity();
  }

  bool MX::isZero() const {
    if (size()==0) {
      return true;
    } else {
      return (*this)->isZero();
    }
  }

  bool MX::isOne() const {
    return (*this)->isOne();
  }

  bool MX::isMinusOne() const {
    return (*this)->isValue(-1);
  }

  bool MX::isTranspose() const {
    return getOp()==OP_TRANSPOSE;
  }

  bool MX::isRegular() const {
    if (isConstant()) {
      return getMatrixValue().isRegular();
    } else {
      casadi_error("Cannot check regularity for symbolic MX");
    }
  }

  MX MX::T() const {
    // Quick return if scalar
    if (isScalar()) {
      return *this;
    } else {
      return (*this)->getTranspose();
    }
  }

  void MX::addToSum(const MX& x) {
    if (isEmpty(true)) {
      *this = x;
    } else {
      *this += x;
    }
  }

  bool MX::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const MXNode*>(ptr)!=0;
  }

  // Helper function
  bool has_empty(const vector<MX>& x, bool both=false) {
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (i->isEmpty(both)) return true;
    }
    return false;
  }

  vector<MX> trim_empty(const vector<MX>& x, bool both=false) {
    vector<MX> ret;
    for (vector<MX>::const_iterator i=x.begin(); i!=x.end(); ++i) {
      if (!i->isEmpty(both)) ret.push_back(*i);
    }
    return ret;
  }

  MX MX::zz_horzcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size1()==ne[0].size1(),
                      "horzcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dimString() <<
                      " and x[0]: " << ne[0].dimString() << ".");
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

  MX MX::zz_diagcat(const vector<MX>& x) {
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
      return x.front()->getDiagcat(x);
    }
  }

  MX MX::zz_vertcat(const vector<MX>& x) {
    // Check dimensions
    if (x.size()>1) {
      vector<MX> ne = trim_empty(x, true);
      for (int i=0;i<ne.size();i++) {
        casadi_assert_message(ne[i].size2()==ne[0].size2(),
                      "vertcat dimension mismatch  " <<
                      "x[" << i << "]:" << ne[i].dimString() <<
                      " and x[0]: " << ne[0].dimString() << ".");
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
    } else if (!x.front().isVector()) {
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

  std::vector<MX> MX::zz_horzsplit(int incr) const {
    casadi_assert(incr>=1);
    vector<int> offset2 = range(0, size2(), incr);
    offset2.push_back(size2());
    return horzsplit(*this, offset2);
  }

  std::vector<MX> MX::zz_diagsplitNative(const std::vector<int>& offset1,
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

    return (*this)->getDiagsplit(offset1, offset2);
  }

  std::vector<MX> MX::zz_vertsplit(const std::vector<int>& offset) const {
    if (isVector()) {
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

  std::vector<MX> MX::zz_vertsplit(int incr) const {
    casadi_assert(incr>=1);
    vector<int> offset1 = range(0, size1(), incr);
    offset1.push_back(size1());
    return vertsplit(*this, offset1);
  }

  std::vector< std::vector<MX> > MX::zz_blocksplit(const std::vector<int>& vert_offset,
                                                   const std::vector<int>& horz_offset) const {
    std::vector<MX> rows = vertsplit(*this, vert_offset);
    std::vector< std::vector<MX> > ret;
    for (int i=0; i<rows.size(); ++i) {
      ret.push_back(horzsplit(rows[i], horz_offset));
    }
    return ret;
  }

  std::vector< std::vector<MX > > MX::zz_blocksplit(int vert_incr, int horz_incr) const {
    casadi_assert(horz_incr>=1);
    casadi_assert(vert_incr>=1);
    vector<int> offset1 = range(0, size1(), vert_incr);
    offset1.push_back(size1());
    vector<int> offset2 = range(0, size2(), horz_incr);
    offset2.push_back(size2());
    return blocksplit(*this, offset1, offset2);
  }

  MX MX::zz_veccat(const vector<MX>& comp) {
    vector<MX> ret = comp;
    for (vector<MX>::iterator i=ret.begin(); i!=ret.end(); ++i) {
      *i = vec(*i);
    }
    return vertcat(ret);
  }

  MX MX::zz_vecNZcat(const vector<MX>& comp) {
    vector<MX> ret = comp;
    for (vector<MX>::iterator i=ret.begin(); i!=ret.end(); ++i) {
      *i = vecNZ(*i);
    }
    return vertcat(ret);
  }

  MX MX::zz_blockcat(const std::vector< std::vector<MX > > &v) {
    // Quick return if no block rows
    if (v.empty()) return MX::sparse(0, 0);

    // Make sure same number of block columns
    int ncols = v.front().size();
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      casadi_assert_message(it->size()==ncols, "blockcat: Inconsistent number of blocl columns");
    }

    // Quick return if no block columns
    if (v.front().empty()) return MX::sparse(0, 0);

    // Horizontally concatenate all columns for each row, then vertically concatenate rows
    std::vector<MX> rows;
    for (vector<vector<MX> >::const_iterator it=v.begin(); it!=v.end(); ++it) {
      rows.push_back(horzcat(*it));
    }
    return vertcat(rows);
  }

  MX MX::zz_blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
    return vertcat(horzcat(A, B), horzcat(C, D));
  }

  MX MX::zz_norm_2() const {
    if (isVector()) {
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

  void MX::zz_simplify() {
    if (!isEmpty(true)) {
      (*this)->simplifyMe(*this);
    }
  }

  MX MX::zz_reshape(std::pair<int, int> rc) const {
    return reshape(*this, rc.first, rc.second);
  }

  MX MX::zz_reshape(int nrow, int ncol) const {
    if (nrow==size1() && ncol==size2())
      return *this;
    else
      return reshape(*this, sparsity().reshape(nrow, ncol));
  }

  MX MX::zz_reshape(const Sparsity& sp) const {
    // quick return if already the right shape
    if (sp==sparsity())
      return *this;

    // make sure that the patterns match
    casadi_assert(sp.isReshape(sparsity()));

    // Create a reshape node
    return (*this)->getReshape(sp);
  }

  MX MX::zz_vec() const {
    if (isVector()) {
      return *this;
    } else {
      return reshape(*this, numel(), 1);
    }
  }

  MX MX::zz_vecNZ() const {
    if (isDense()) {
      return vec(*this);
    } else {
      return (*this)->getGetNonzeros(Sparsity::dense(size(), 1), range(size()));
    }
  }

  MX MX::zz_if_else(const MX &if_true, const MX &if_false) const {
    return casadi::if_else_zero(*this, if_true) + casadi::if_else_zero(!*this, if_false);
  }

  MX MX::zz_unite(const MX& B) const {
    // Join the sparsity patterns
    std::vector<unsigned char> mapping;
    Sparsity sp = sparsity().patternUnion(B.sparsity(), mapping);

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

  MX MX::zz_repmat(int n, int m) const {
    // Quick return if possible
    if (n==1 &&  m==1) return *this;

    // First concatenate horizontally
    MX col = horzcat(std::vector<MX >(m, *this));

    // Then vertically
    return vertcat(std::vector<MX >(n, col));
  }

  MX MX::zz_dense() const {
    MX ret = *this;
    ret.densify();
    return ret;
  }

  MX MX::zz_createParent(std::vector<MX> &deps) {
    // First check if arguments are symbolic
    for (int k=0;k<deps.size();k++) {
      if (!deps[k].isSymbolic())
          throw CasadiException("createParent: the argumenst must be pure symbolic");
    }

    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1, 0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }

    // Create the parent
    MX P = MX::sym("P", index[deps.size()], 1);

    std::vector<MX> Ps = vertsplit(P, index);

    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      deps[k] = MX(deps[k].sparsity(), Ps[k]);
    }

    return P;
  }

  MX MX::zz_createParent(const std::vector<Sparsity> &deps, std::vector<MX>& children) {
    // Collect the sizes of the depenencies
    std::vector<int> index(deps.size()+1, 0);
    for (int k=0;k<deps.size();k++) {
      index[k+1] =  index[k] + deps[k].size();
    }

    // Create the parent
    MX P = MX::sym("P", index[deps.size()], 1);
    std::vector<MX> Ps = vertsplit(P, index);

    children.resize(deps.size());

    // Make the arguments dependent on the parent
    for (int k=0;k<deps.size();k++) {
      children[k] =  MX(deps[k], Ps[k]);
    }

    return P;
  }

  MX MX::zz_createParent(const std::vector<MX> &deps, std::vector<MX>& children) {
    children = deps;
    MX P = createParent(children);
    return P;
  }

  MX MX::zz_diag() const {
    // Nonzero mapping
    std::vector<int> mapping;

    // Get the sparsity
    Sparsity sp = sparsity().getDiag(mapping);

    // Create a reference to the nonzeros
    return (*this)->getGetNonzeros(sp, mapping);
  }

  MX MX::zz_blkdiag(const std::vector<MX> &A) {
    return diagcat(A);
  }

  MX MX::zz_blkdiag(const MX& B) const {
    std::vector<MX> ret;
    ret.push_back(*this);
    ret.push_back(B);
    return blkdiag(ret);
  }

  int MX::zz_countNodes() const {
    MXFunction f(vector<MX>(), *this);
    f.init();
    return f.countNodes();
  }

  MX MX::zz_sumCols() const {
    return mul(*this, MX::ones(size2(), 1));
  }

  MX MX::zz_sumRows() const {
    return mul(MX::ones(1, size1()), *this);
  }

  MX MX::zz_sumAll() const {
    return sumRows(sumCols(*this));
  }


  MX MX::zz_polyval(const MX& x) const {
    casadi_assert_message(isDense(), "polynomial coefficients vector must be a vector");
    casadi_assert_message(isVector() && size()>0, "polynomial coefficients must be a vector");
    MX ret = (*this)[0];
    for (int i=1; i<size(); ++i) {
      ret = ret*x + (*this)[i];
    }
    return ret;
  }

  std::string MX::zz_getOperatorRepresentation(const std::vector<std::string>& args) const {
    std::stringstream s;
    const MXNode* node = dynamic_cast<const MXNode*>(get());
    node->printPart(s, 0);
    if (isUnary()) {
      s << args[0];
      node->printPart(s, 1);
    } else {
      for (int i=0;i<args.size();++i) {
        s << args[i];
        node->printPart(s, 1+i);
      }
    }

    return s.str();
  }

  void MX::zz_substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, bool reverse) {
    // Empty vector
    vector<MX> ex;
    substituteInPlace(v, vdef, ex, reverse);
  }

  void MX::zz_substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef,
                         std::vector<MX>& ex, bool reverse) {
    casadi_assert_message(v.size()==vdef.size(),
                          "Mismatch in the number of expression to substitute.");
    for (int k=0; k<v.size(); ++k) {
      casadi_assert_message(v[k].isSymbolic(), "Variable " << k << " is not symbolic");
      casadi_assert_message(v[k].sparsity() == vdef[k].sparsity(),
                            "Inconsistent sparsity for variable " << k << ".");
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
    MXFunction f(f_in, f_out);
    f.init();

    // Get references to the internal data structures
    std::vector<MXAlgEl>& algorithm = f->algorithm_;
    vector<MX> work(f.getWorkSize());
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

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
          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i];
            input_p[i] = el<0 ? 0 : &work.at(el);
          }

          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i];
            output_p[i] = el<0 ? 0 : &work.at(el);
          }

          it->data->evaluateMX(input_p, output_p, dummy_p, dummy_p, dummy_p, dummy_p, false);
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
      if (!v[k].isEqual(vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Otherwise, evaluate symbolically
    MXFunction F(v, ex);
    F.init();
    return F.call(vdef, true);
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
    MXFunction f(vector<MX>(), ex);
    f.init();

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = f.algorithm();
    vector<MX> swork(f.getWorkSize());

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
    vector<MX> f_out(f.getNumOutputs());

    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

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

          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i];
            if (el>=0) node_tainted =  node_tainted || tainted[el];
            input_p[i] = el<0 ? 0 : &swork[el];
          }

          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i];
            output_p[i] = el<0 ? 0 : &swork[el];
            if (el>=0) tainted[el] = node_tainted;
          }

          if (it->res.size()==1 && it->res[0]>=0 && !node_tainted) {
            int el = it->res[0];
            swork[el] = it->data;
          } else {
            const_cast<MX&>(it->data)->evaluateMX(input_p, output_p,
                                                  dummy_p, dummy_p, dummy_p, dummy_p, false);
          }
        }
      }
    }

    bool all_found=true;
    for (int i=0;i<expr.size();++i) {
      all_found = all_found && expr_found[i];
    }

    //casadi_assert_message(all_found,
    //             "MXFunctionInternal::extractNodes(const std::vector<MX>& expr):"
    //             " failed to locate all input expr."
    //             << std::endl << "Here's a boolean list showing which ones where found: "
    //             << expr_found);

    return f_out;

  }

  void MX::zz_extractShared(std::vector<MX>& ex, std::vector<MX>& v, std::vector<MX>& vdef,
                            const std::string& v_prefix, const std::string& v_suffix) {

    // Sort the expression
    MXFunction f(vector<MX>(), ex);
    f.init();

    // Get references to the internal data structures
    const vector<MXAlgEl>& algorithm = f.algorithm();
    vector<MX> work(f.getWorkSize());

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
    MXPtrV input_p, output_p;
    MXPtrVV dummy_p;

    // Evaluate the algorithm
    k=0;
    for (vector<MXAlgEl>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it, ++k) {
      switch (it->op) {
      case OP_OUTPUT:     ex[it->res.front()] = work[it->arg.front()];      break;
      case OP_CONST:
      case OP_PARAMETER:  work[it->res.front()] = it->data; break;
      default:
        {
          // Pointers to the arguments of the evaluation
          input_p.resize(it->arg.size());
          for (int i=0; i<input_p.size(); ++i) {
            int el = it->arg[i]; // index of the argument
            input_p[i] = el<0 ? 0 : &work[el];
          }

          // Pointers to the result of the evaluation
          output_p.resize(it->res.size());
          for (int i=0; i<output_p.size(); ++i) {
            int el = it->res[i]; // index of the output
            output_p[i] = el<0 ? 0 : &work[el];
          }

          // Evaluate atomic operation
          const_cast<MX&>(it->data)->evaluateMX(input_p, output_p,
                                                dummy_p, dummy_p, dummy_p, dummy_p, false);

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

  void MX::zz_printCompact(std::ostream &stream) const {
    // Extract shared subexpressions from ex
    vector<MX> v, vdef;
    vector<MX> ex_extracted(1, *this);
    extractShared(ex_extracted, v, vdef, "@", "");

    // Print the expression without shared subexpressions
    ex_extracted.front().print(stream);

    // Print the shared subexpressions
    if (!v.empty()) {
      stream << endl << "where:" << endl;
      for (int i=0; i<v.size(); ++i) {
        stream << v[i] << " := " << vdef[i] << endl;
      }
    }
  }

  MX MX::zz_jacobian(const MX &arg) const {
    MXFunction temp(arg, *this); // make a runtime
    temp.setOption("name", "helper_jacobian_MX");
    temp.init();
    return temp.jac();
  }

  MX MX::zz_gradient(const MX &arg) const {
    MXFunction temp(arg, *this); // make a runtime
    temp.setOption("name", "helper_gradient_MX");
    temp.init();
    return temp.grad();
  }

  MX MX::zz_tangent(const MX &arg) const {
    MXFunction temp(arg, *this); // make a runtime
    temp.setOption("name", "helper_tangent_MX");
    temp.init();
    return temp.tang();
  }

  MX MX::zz_det() const {
    return (*this)->getDeterminant();
  }

  MX MX::zz_inv() const {
    return (*this)->getInverse();
  }

  std::vector<MX> MX::zz_getSymbols() const {
    MXFunction f(std::vector<MX>(), *this);
    f.init();
    return f.getFree();
  }

  std::vector<MX> MX::zz_getSymbols(const std::vector<MX>& e) {
    MXFunction f(std::vector<MX>(), e);
    f.init();
    return f.getFree();
  }

  bool MX::zz_dependsOn(const std::vector<MX> &arg) const {
    if (size()==0) return false;

    // Construct a temporary algorithm
    MXFunction temp(arg, *this);
    temp.init();
    temp.spInit(true);

    for (int i=0;i<temp.getNumInputs();++i) {
      bvec_t* input_ =  get_bvec_t(temp.input(i).data());
      std::fill(input_, input_+temp.input(i).size(), bvec_t(1));
    }
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.output().size();++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  MX MX::zz_matrix_expand(const MX& e, const std::vector<MX> &boundary) {
    std::vector<MX> e_v(1, e);
    return matrix_expand(e_v, boundary).at(0);
  }

  std::vector<MX> MX::zz_matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary) {

    // Create symbols for boundary nodes
    std::vector<MX> syms(boundary.size());

    for (int i=0;i<syms.size();++i) {
      syms[i] = MX::sym("x", boundary[i].sparsity());
    }

    // Substitute symbols for boundary nodes
    std::vector<MX> ret = graph_substitute(e, boundary, syms);

    // Obtain list of dependents
    std::vector<MX> v = getSymbols(ret);

    // Construct an MXFunction with it
    MXFunction f(v, ret);
    f.init();

    // Expand to SXFunction
    SXFunction s = f.expand();
    s.init();

    return s.call(graph_substitute(v, syms, boundary), true);
  }

  MX MX::zz_kron(const MX& b) const {
    const Sparsity &a_sp = sparsity();
    MX filler = MX::sparse(b.shape());
    std::vector< std::vector< MX > > blocks(size1(), std::vector< MX >(size2(), filler));
    for (int i=0; i<size1(); ++i) {
      for (int j=0; j<size2(); ++j) {
        int k = a_sp.elem(i, j);
        if (k!=-1) {
          blocks[i][j] = (*this)[k]*b;
        }
      }
    }
    return blockcat(blocks);
  }

  MX MX::zz_solve(const MX& b, const std::string& lsolver, const Dictionary& dict) const {
    LinearSolver mysolver(lsolver, sparsity(), b.size2());
    mysolver.setOption(dict);
    mysolver.init();
    return mysolver.solve(*this, b, false);
  }

  MX MX::zz_pinv(const std::string& lsolver, const Dictionary& dict) const {
    if (size1()>=size2()) {
      return solve(mul(T(), *this), T(), lsolver, dict);
    } else {
      return solve(mul(*this, T()), *this, lsolver, dict).T();
    }
  }

  MX MX::zz_nullspace() const {
    SX n = SX::sym("A", sparsity());
    SXFunction f(n, nullspace(n));
    f.setOption("name", "nullspace");
    f.init();
    return f(*this).at(0);
  }

} // namespace casadi
