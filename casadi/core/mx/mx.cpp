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

  const MX MX::sub(int j, const std::vector<int>& i) const {
    return sub(std::vector<int>(1, j), i);
  }

  const MX MX::sub(const std::vector<int>& jj, const std::vector<int>& ii) const {
    // Nonzero mapping from submatrix to full
    std::vector<int> mapping;

    // Get the sparsity pattern
    Sparsity sp = sparsity().sub(jj, ii, mapping);

    // Create return MX
    return (*this)->getGetNonzeros(sp, mapping);
  }

  const MX MX::sub(int j, int i) const {
    int ind = sparsity().getNZ(j, i);
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

  MX MX::mul_full(const MX& y, const Sparsity &z) const {
    const MX& x = *this;
    return x->getMultiplication(y, z);
  }

  MX MX::mul(const MX& y, const Sparsity &z) const {
    return mul_smart(y, z);
  }

  MX MX::inner_prod(const MX& y) const {
    return (*this)->getInnerProd(y);
  }

  MX MX::outer_prod(const MX& y) const {
    return mul(y.T());
  }

  MX MX::__pow__(const MX& n) const {
    if (n->getOp()==OP_CONST) {
      return MX::binary(OP_CONSTPOW, *this, n);
    } else {
      return MX::binary(OP_POW, *this, n);
    }
  }

  MX MX::constpow(const MX& b) const {
    return binary(OP_CONSTPOW, *this, b);
  }

  MX MX::fmin(const MX& b) const {
    return binary(OP_FMIN, *this, b);
  }

  MX MX::fmax(const MX& b) const {
    return binary(OP_FMAX, *this, b);
  }

  MX MX::fmod(const MX& b) const {
    return binary(OP_FMOD, *this, b);
  }

  MX MX::arctan2(const MX& b) const {
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

  MX MX::exp() const {
    return (*this)->getUnary(OP_EXP);
  }

  MX MX::log() const {
    return (*this)->getUnary(OP_LOG);
  }

  MX MX::log10() const {
    return log()*(1/std::log(10.));
  }

  MX MX::sqrt() const {
    return (*this)->getUnary(OP_SQRT);
  }

  MX MX::sin() const {
    return (*this)->getUnary(OP_SIN);
  }

  MX MX::cos() const {
    return (*this)->getUnary(OP_COS);
  }

  MX MX::tan() const {
    return (*this)->getUnary(OP_TAN);
  }

  MX MX::arcsin() const {
    return (*this)->getUnary(OP_ASIN);
  }

  MX MX::arccos() const {
    return (*this)->getUnary(OP_ACOS);
  }

  MX MX::arctan() const {
    return (*this)->getUnary(OP_ATAN);
  }

  MX MX::sinh() const {
    return (*this)->getUnary(OP_SINH);
  }

  MX MX::cosh() const {
    return (*this)->getUnary(OP_COSH);
  }

  MX MX::tanh() const {
    return (*this)->getUnary(OP_TANH);
  }

  MX MX::arcsinh() const {
    return (*this)->getUnary(OP_ASINH);
  }

  MX MX::arccosh() const {
    return (*this)->getUnary(OP_ACOSH);
  }

  MX MX::arctanh() const {
    return (*this)->getUnary(OP_ATANH);
  }

  MX MX::floor() const {
    return (*this)->getUnary(OP_FLOOR);
  }

  MX MX::ceil() const {
    return (*this)->getUnary(OP_CEIL);
  }

  MX MX::fabs() const {
    return (*this)->getUnary(OP_FABS);
  }

  MX MX::sign() const {
    return (*this)->getUnary(OP_SIGN);
  }

  MX MX::__copysign__(const MX& y) const {
    return MX::binary(OP_COPYSIGN, *this, y);
  }

  MX MX::erfinv() const {
    return (*this)->getUnary(OP_ERFINV);
  }

  MX MX::erf() const {
    return (*this)->getUnary(OP_ERF);
  }

  MX MX::logic_not() const {
    return (*this)->getUnary(OP_NOT);
  }

  void MX::lift(const MX& x_guess) {
    casadi_assert(sparsity()==x_guess.sparsity());
    *this = (*this)->getBinary(OP_LIFT, x_guess, false, false);
  }

  MX MX::__add__(const MX& y) const {
    return MX::binary(OP_ADD, *this, y);
  }

  MX MX::__sub__(const MX& y) const {
    return MX::binary(OP_SUB, *this, y);
  }

  MX MX::__mul__(const MX& y) const {
    return MX::binary(OP_MUL, *this, y);
  }

  MX MX::__div__(const MX& y) const {
    return MX::binary(OP_DIV, *this, y);
  }

  MX MX::__lt__(const MX& y) const {
    return MX::binary(OP_LT, *this, y);
  }

  MX MX::__le__(const MX& y) const {
    return MX::binary(OP_LE, *this, y);
  }

  MX MX::__eq__(const MX& y) const {
    return MX::binary(OP_EQ, *this, y);
  }

  MX MX::__ne__(const MX& y) const {
    return MX::binary(OP_NE, *this, y);
  }

  MX MX::logic_and(const MX& y) const {
    return MX::binary(OP_AND, *this, y);
  }

  MX MX::logic_or(const MX& y) const {
    return MX::binary(OP_OR, *this, y);
  }

  MX MX::if_else_zero(const MX& y) const {
    return MX::binary(OP_IF_ELSE_ZERO, *this, y);
  }

  MX MX::__constpow__(const MX& b) const { return (*this).constpow(b);}
  MX MX::__mrdivide__(const MX& b) const
  { if (b.isScalar()) return *this/b; throw CasadiException("mrdivide: Not implemented");}
  MX MX::__mpower__(const MX& b) const
  { return pow(*this, b); throw CasadiException("mpower: Not implemented");}

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

  MX MX::trans() const {
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

} // namespace casadi
