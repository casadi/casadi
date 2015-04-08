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


#include "sx_element.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/generic_expression_tools.hpp"
#include <stack>
#include <cassert>
#include "../casadi_math.hpp"
#include "constant_sx.hpp"
#include "symbolic_sx.hpp"
#include "unary_sx.hpp"
#include "binary_sx.hpp"
#include "../casadi_options.hpp"
#include "../function/sx_function_internal.hpp"
#include "sx_tools.hpp"

using namespace std;
namespace casadi {

  // Allocate storage for the caching
  CACHING_MAP<int, IntegerSX*> IntegerSX::cached_constants_;
  CACHING_MAP<double, RealtypeSX*> RealtypeSX::cached_constants_;

  SXElement::SXElement() {
    node = casadi_limits<SXElement>::nan.node;
    node->count++;
  }

  SXElement::SXElement(SXNode* node_, bool dummy) : node(node_) {
    node->count++;
  }

  SXElement SXElement::create(SXNode* node) {
    return SXElement(node, false);
  }

  SXElement::SXElement(const SXElement& scalar) {
    node = scalar.node;
    node->count++;
  }

  SXElement::SXElement(double val) {
    int intval = static_cast<int>(val);
    if (val-intval == 0) { // check if integer
      if (intval == 0)             node = casadi_limits<SXElement>::zero.node;
      else if (intval == 1)        node = casadi_limits<SXElement>::one.node;
      else if (intval == 2)        node = casadi_limits<SXElement>::two.node;
      else if (intval == -1)       node = casadi_limits<SXElement>::minus_one.node;
      else                        node = IntegerSX::create(intval);
      node->count++;
    } else {
      if (isnan(val))              node = casadi_limits<SXElement>::nan.node;
      else if (isinf(val))         node = val > 0 ? casadi_limits<SXElement>::inf.node :
                                      casadi_limits<SXElement>::minus_inf.node;
      else                        node = RealtypeSX::create(val);
      node->count++;
    }
  }

  SXElement SXElement::sym(const std::string& name) {
    return create(new SymbolicSX(name));
  }

  SXElement::~SXElement() {
    if (--node->count == 0) delete node;
  }

  SXElement& SXElement::operator=(const SXElement &scalar) {
    // quick return if the old and new pointers point to the same object
    if (node == scalar.node) return *this;

    // decrease the counter and delete if this was the last pointer
    if (--node->count == 0) delete node;

    // save the new pointer
    node = scalar.node;
    node->count++;
    return *this;
  }

  void SXElement::assignIfDuplicate(const SXElement& scalar, int depth) {
    casadi_assert(depth>=1);
    if (!isEqual(*this, scalar, 0) && isEqual(*this, scalar, depth)) {
      *this = scalar;
    }
  }

  SXNode* SXElement::assignNoDelete(const SXElement& scalar) {
    // Return value
    SXNode* ret = node;

    // quick return if the old and new pointers point to the same object
    if (node == scalar.node) return ret;

    // decrease the counter but do not delete if this was the last pointer
    --node->count;

    // save the new pointer
    node = scalar.node;
    node->count++;

    // Return a pointer to the old node
    return ret;
  }

  SXElement& SXElement::operator=(double scalar) {
    return *this = SXElement(scalar);
  }

  void SXElement::repr(std::ostream &stream, bool trailing_newline) const {
    print(stream, trailing_newline);
  }

  void SXElement::print(std::ostream &stream, bool trailing_newline) const {
    node->print(stream);
    if (trailing_newline) stream << std::endl;
  }

  void SXElement::print(std::ostream &stream, long& remaining_calls) const {
    if (remaining_calls>0) {
      remaining_calls--;
      node->print(stream, remaining_calls);
    } else {
      stream << "...";
    }
  }

  SXElement SXElement::operator-() const {
    if (isOp(OP_NEG))
      return getDep();
    else if (isZero())
      return 0;
    else if (isMinusOne())
      return 1;
    else if (isOne())
      return -1;
    else
      return UnarySX::create(OP_NEG, *this);
  }

  SXElement SXElement::zz_sign() const {
    return UnarySX::create(OP_SIGN, *this);
  }

  SXElement SXElement::__copysign__(const SXElement &y) const {
    return BinarySX::create(OP_COPYSIGN, *this, y);
  }

  SXElement SXElement::zz_erfinv() const {
    return UnarySX::create(OP_ERFINV, *this);
  }

  bool SXElement::__nonzero__() const {
    if (isConstant()) return !isZero();
    casadi_error("Cannot compute the truth value of a CasADi SXElement symbolic expression.")
  }

  template<>
  bool Matrix<SXElement>::__nonzero__() const {
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you "
                                  "provided a shape" << dimString());}
    return at(0).__nonzero__();
  }

  SXElement SXElement::zz_plus(const SXElement& y) const {
    // NOTE: Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_ADD, *this, y);

    if (isZero())
      return y;
    else if (y->isZero()) // term2 is zero
      return *this;
    else if (y.isOp(OP_NEG))  // x + (-y) -> x - y
      return zz_minus(-y);
    else if (isOp(OP_NEG)) // (-x) + y -> y - x
      return y.zz_minus(getDep());
    else if (isOp(OP_MUL) && y.isOp(OP_MUL) &&
            getDep(0).isConstant() && getDep(0).getValue()==0.5 &&
            y.getDep(0).isConstant() && y.getDep(0).getValue()==0.5 &&
             isEqual(y.getDep(1), getDep(1), SXNode::eq_depth_)) // 0.5x+0.5x = x
      return getDep(1);
    else if (isOp(OP_DIV) && y.isOp(OP_DIV) &&
             getDep(1).isConstant() && getDep(1).getValue()==2 &&
             y.getDep(1).isConstant() && y.getDep(1).getValue()==2 &&
             isEqual(y.getDep(0), getDep(0), SXNode::eq_depth_)) // x/2+x/2 = x
      return getDep(0);
    else if (isOp(OP_SUB) && isEqual(getDep(1), y, SXNode::eq_depth_))
      return getDep(0);
    else if (y.isOp(OP_SUB) && isEqual(*this, y.getDep(1), SXNode::eq_depth_))
      return y.getDep(0);
    else if (isOp(OP_SQ) && y.isOp(OP_SQ) &&
             ((getDep().isOp(OP_SIN) && y.getDep().isOp(OP_COS))
              || (getDep().isOp(OP_COS) && y.getDep().isOp(OP_SIN)))
             && isEqual(getDep().getDep(), y.getDep().getDep(), SXNode::eq_depth_))
      return 1; // sin^2 + cos^2 -> 1
    else // create a new branch
      return BinarySX::create(OP_ADD, *this, y);
  }

  SXElement SXElement::zz_minus(const SXElement& y) const {
    // Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_SUB, *this, y);

    if (y->isZero()) // term2 is zero
      return *this;
    if (isZero()) // term1 is zero
      return -y;
    if (isEqual(*this, y, SXNode::eq_depth_)) // the terms are equal
      return 0;
    else if (y.isOp(OP_NEG)) // x - (-y) -> x + y
      return *this + y.getDep();
    else if (isOp(OP_ADD) && isEqual(getDep(1), y, SXNode::eq_depth_))
      return getDep(0);
    else if (isOp(OP_ADD) && isEqual(getDep(0), y, SXNode::eq_depth_))
      return getDep(1);
    else if (y.isOp(OP_ADD) && isEqual(*this, y.getDep(1), SXNode::eq_depth_))
      return -y.getDep(0);
    else if (y.isOp(OP_ADD) && isEqual(*this, y.getDep(0), SXNode::eq_depth_))
      return -y.getDep(1);
    else if (isOp(OP_NEG))
      return -(getDep() + y);
    else // create a new branch
      return BinarySX::create(OP_SUB, *this, y);
  }

  SXElement SXElement::zz_times(const SXElement& y) const {

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_MUL, *this, y);

    // Only simplifications that do not result in extra nodes area allowed
    if (isEqual(y, *this, SXNode::eq_depth_))
      return sq();
    else if (!isConstant() && y.isConstant())
      return y.zz_times(*this);
    else if (isZero() || y->isZero()) // one of the terms is zero
      return 0;
    else if (isOne()) // term1 is one
      return y;
    else if (y->isOne()) // term2 is one
      return *this;
    else if (y->isMinusOne())
      return -(*this);
    else if (isMinusOne())
      return -y;
    else if (y.isOp(OP_INV))
      return (*this)/y.inv();
    else if (isOp(OP_INV))
      return y/inv();
    else if (isConstant() && y.isOp(OP_MUL) && y.getDep(0).isConstant() &&
            getValue()*y.getDep(0).getValue()==1) // 5*(0.2*x) = x
      return y.getDep(1);
    else if (isConstant() && y.isOp(OP_DIV) && y.getDep(1).isConstant() &&
            getValue()==y.getDep(1).getValue()) // 5*(x/5) = x
      return y.getDep(0);
    else if (isOp(OP_DIV) && isEqual(getDep(1), y, SXNode::eq_depth_)) // ((2/x)*x)
      return getDep(0);
    else if (y.isOp(OP_DIV) &&
             isEqual(y.getDep(1), *this, SXNode::eq_depth_)) // ((2/x)*x)
      return y.getDep(0);
    else if (isOp(OP_NEG))
      return -(getDep() * y);
    else if (y.isOp(OP_NEG))
      return -(*this * y.getDep());
    else     // create a new branch
      return BinarySX::create(OP_MUL, *this, y);
  }


  bool SXElement::isDoubled() const {
    return isOp(OP_ADD) && isEqual(getDep(0), getDep(1), SXNode::eq_depth_);
  }

  SXElement SXElement::zz_rdivide(const SXElement& y) const {
    // Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_DIV, *this, y);

    if (y->isZero()) // term2 is zero
      return casadi_limits<SXElement>::nan;
    else if (isZero()) // term1 is zero
      return 0;
    else if (y->isOne()) // term2 is one
      return *this;
    else if (y->isMinusOne())
      return -(*this);
    else if (isEqual(*this, y, SXNode::eq_depth_)) // terms are equal
      return 1;
    else if (isDoubled() && isEqual(y, 2))
      return getDep(0);
    else if (isOp(OP_MUL) && isEqual(y, getDep(0), SXNode::eq_depth_))
      return getDep(1);
    else if (isOp(OP_MUL) && isEqual(y, getDep(1), SXNode::eq_depth_))
      return getDep(0);
    else if (isOne())
      return y.inv();
    else if (y.isOp(OP_INV))
      return (*this)*y.inv();
    else if (isDoubled() && y.isDoubled())
      return getDep(0) / y->dep(0);
    else if (y.isConstant() && isOp(OP_DIV) && getDep(1).isConstant() &&
            y.getValue()*getDep(1).getValue()==1) // (x/5)/0.2
      return getDep(0);
    else if (y.isOp(OP_MUL) &&
             isEqual(y.getDep(1), *this, SXNode::eq_depth_)) // x/(2*x) = 1/2
      return BinarySX::create(OP_DIV, 1, y.getDep(0));
    else if (isOp(OP_NEG) &&
             isEqual(getDep(0), y, SXNode::eq_depth_))      // (-x)/x = -1
      return -1;
    else if (y.isOp(OP_NEG) &&
             isEqual(y.getDep(0), *this, SXNode::eq_depth_))      // x/(-x) = 1
      return -1;
    else if (y.isOp(OP_NEG) && isOp(OP_NEG) &&
             isEqual(getDep(0), y.getDep(0), SXNode::eq_depth_))  // (-x)/(-x) = 1
      return 1;
    else if (isOp(OP_DIV) && isEqual(y, getDep(0), SXNode::eq_depth_))
      return getDep(1).inv();
    else if (isOp(OP_NEG))
      return -(getDep() / y);
    else if (y.isOp(OP_NEG))
      return -(*this / y.getDep());
    else // create a new branch
      return BinarySX::create(OP_DIV, *this, y);
  }

  SXElement SXElement::inv() const {
    if (isOp(OP_INV)) {
      return getDep(0);
    } else {
      return UnarySX::create(OP_INV, *this);
    }
  }

  SX SXElement::zz_min(const SX& b) const {
    return fmin(SX(*this), b);
  }
  SX SXElement::zz_max(const SX& b) const {
    return fmax(SX(*this), b);
  }
  SX SXElement::constpow(const SX& n) const {
    return SX(*this).__constpow__(n);
  }
  SX SXElement::__copysign__(const SX& n) const {
    return SX(*this).__copysign__(n);
  }

  SX SXElement::zz_atan2(const SX& b) const {
    return atan2(SX(*this), b);
  }

  SXElement SXElement::zz_le(const SXElement& y) const {
    if ((y-(*this)).isNonNegative())
      return 1;
    else
      return BinarySX::create(OP_LE, *this, y);
  }

  SXElement SXElement::zz_lt(const SXElement& y) const {
    if (((*this)-y).isNonNegative())
      return 0;
    else
      return BinarySX::create(OP_LT, *this, y);
  }

  SXElement SXElement::zz_eq(const SXElement& y) const {
    if (isEqual(*this, y))
      return 1;
    else
      return BinarySX::create(OP_EQ, *this, y);
  }

  SXElement SXElement::zz_ne(const SXElement& y) const {
    if (isEqual(*this, y))
      return 0;
    else
      return BinarySX::create(OP_NE, *this, y);
  }

  SXNode* SXElement::get() const {
    return node;
  }

  const SXNode* SXElement::operator->() const {
    return node;
  }

  SXNode* SXElement::operator->() {
    return node;
  }

  SXElement if_else(const SXElement& cond, const SXElement& if_true, const SXElement& if_false) {
    return if_else_zero(cond, if_true) + if_else_zero(!cond, if_false);
  }

  SXElement SXElement::binary(int op, const SXElement& x, const SXElement& y) {
    return BinarySX::create(Operation(op), x, y);
  }

  SXElement SXElement::unary(int op, const SXElement& x) {
    return UnarySX::create(Operation(op), x);
  }

  bool SXElement::isLeaf() const {
    if (!node) return true;
    return isConstant() || isSymbolic();
  }

  bool SXElement::isCommutative() const {
    if (!hasDep()) throw CasadiException("SX::isCommutative: must be binary");
    return operation_checker<CommChecker>(getOp());
  }

  bool SXElement::isConstant() const {
    return node->isConstant();
  }

  bool SXElement::isInteger() const {
    return node->isInteger();
  }

  bool SXElement::isSymbolic() const {
    return node->isSymbolic();
  }

  bool SXElement::hasDep() const {
    return node->hasDep();
  }

  bool SXElement::isZero() const {
    return node->isZero();
  }

  bool SXElement::isAlmostZero(double tol) const {
    return node->isAlmostZero(tol);
  }

  bool SXElement::isOne() const {
    return node->isOne();
  }

  bool SXElement::isMinusOne() const {
    return node->isMinusOne();
  }

  bool SXElement::isNan() const {
    return node->isNan();
  }

  bool SXElement::isInf() const {
    return node->isInf();
  }

  bool SXElement::isMinusInf() const {
    return node->isMinusInf();
  }

  const std::string& SXElement::getName() const {
    return node->getName();
  }

  int SXElement::getOp() const {
    return node->getOp();
  }

  bool SXElement::isOp(int op) const {
    return hasDep() && op==getOp();
  }

  bool SXElement::zz_isEqual(const SXElement& ex, int depth) const {
    if (node==ex.get())
      return true;
    else if (depth>0)
      return node->zz_isEqual(ex.get(), depth);
    else
      return false;
  }

  bool SXElement::isNonNegative() const {
    if (isConstant())
      return getValue()>=0;
    else if (isOp(OP_SQ) || isOp(OP_FABS))
      return true;
    else
      return false;
  }

  double SXElement::getValue() const {
    return node->getValue();
  }

  int SXElement::getIntValue() const {
    return node->getIntValue();
  }

  SXElement SXElement::getDep(int ch) const {
    casadi_assert(ch==0 || ch==1;)
      return node->dep(ch);
  }

  int SXElement::getNdeps() const {
    if (!hasDep()) throw CasadiException("SX::getNdeps: must be binary");
    return casadi_math<double>::ndeps(getOp());
  }

  long SXElement::__hash__() const {
    if (!node) return 0;
    return (long) node;
  }

  // node corresponding to a constant 0
  const SXElement casadi_limits<SXElement>::zero(new ZeroSX(), false);
  // node corresponding to a constant 1
  const SXElement casadi_limits<SXElement>::one(new OneSX(), false);
  // node corresponding to a constant 2
  const SXElement casadi_limits<SXElement>::two(IntegerSX::create(2), false);
  // node corresponding to a constant -1
  const SXElement casadi_limits<SXElement>::minus_one(new MinusOneSX(), false);
  const SXElement casadi_limits<SXElement>::nan(new NanSX(), false);
  const SXElement casadi_limits<SXElement>::inf(new InfSX(), false);
  const SXElement casadi_limits<SXElement>::minus_inf(new MinusInfSX(), false);

  bool casadi_limits<SXElement>::isZero(const SXElement& val) {
    return val.isZero();
  }

  bool casadi_limits<SXElement>::isAlmostZero(const SXElement& val, double tol) {
    return val.isAlmostZero(tol);
  }

  bool casadi_limits<SXElement>::isOne(const SXElement& val) {
    return val.isOne();
  }

  bool casadi_limits<SXElement>::isMinusOne(const SXElement& val) {
    return val.isMinusOne();
  }

  bool casadi_limits<SXElement>::isConstant(const SXElement& val) {
    return val.isConstant();
  }

  bool casadi_limits<SXElement>::isInteger(const SXElement& val) {
    return val.isInteger();
  }

  bool casadi_limits<SXElement>::isInf(const SXElement& val) {
    return val.isInf();
  }

  bool casadi_limits<SXElement>::isMinusInf(const SXElement& val) {
    return val.isMinusInf();
  }

  bool casadi_limits<SXElement>::isNaN(const SXElement& val) {
    return val.isNan();
  }

  SXElement SXElement::zz_exp() const {
    return UnarySX::create(OP_EXP, *this);
  }

  SXElement SXElement::zz_log() const {
    return UnarySX::create(OP_LOG, *this);
  }

  SXElement SXElement::zz_log10() const {
    return log(*this)*(1/std::log(10.));
  }

  SXElement SXElement::zz_sqrt() const {
    if (isOp(OP_SQ))
      return fabs(getDep());
    else
      return UnarySX::create(OP_SQRT, *this);
  }

  SXElement SXElement::sq() const {
    if (isOp(OP_SQRT))
      return getDep();
    else if (isOp(OP_NEG))
      return getDep().sq();
    else
      return UnarySX::create(OP_SQ, *this);
  }

  SXElement SXElement::zz_sin() const {
    return UnarySX::create(OP_SIN, *this);
  }

  SXElement SXElement::zz_cos() const {
    return UnarySX::create(OP_COS, *this);
  }

  SXElement SXElement::zz_tan() const {
    return UnarySX::create(OP_TAN, *this);
  }

  SXElement SXElement::zz_asin() const {
    return UnarySX::create(OP_ASIN, *this);
  }

  SXElement SXElement::zz_acos() const {
    return UnarySX::create(OP_ACOS, *this);
  }

  SXElement SXElement::zz_atan() const {
    return UnarySX::create(OP_ATAN, *this);
  }

  SXElement SXElement::zz_sinh() const {
    if (isZero())
      return 0;
    else
      return UnarySX::create(OP_SINH, *this);
  }

  SXElement SXElement::zz_cosh() const {
    if (isZero())
      return 1;
    else
      return UnarySX::create(OP_COSH, *this);
  }

  SXElement SXElement::zz_tanh() const {
    if (isZero())
      return 0;
    else
      return UnarySX::create(OP_TANH, *this);
  }

  SXElement SXElement::zz_atanh() const {
    if (isZero())
      return 0;
    else
      return UnarySX::create(OP_ATANH, *this);
  }

  SXElement SXElement::zz_acosh() const {
    if (isOne())
      return 0;
    else
      return UnarySX::create(OP_ACOSH, *this);
  }

  SXElement SXElement::zz_asinh() const {
    if (isZero())
      return 0;
    else
      return UnarySX::create(OP_ASINH, *this);
  }

  SXElement SXElement::zz_floor() const {
    return UnarySX::create(OP_FLOOR, *this);
  }

  SXElement SXElement::zz_ceil() const {
    return UnarySX::create(OP_CEIL, *this);
  }

  SXElement SXElement::zz_mod(const SXElement &b) const {
    return BinarySX::create(OP_FMOD, *this, b);
  }

  SXElement SXElement::zz_erf() const {
    return UnarySX::create(OP_ERF, *this);
  }

  SXElement SXElement::zz_abs() const {
    if (isOp(OP_FABS) || isOp(OP_SQ))
      return *this;
    else
      return UnarySX::create(OP_FABS, *this);
  }

  SXElement::operator SX() const {
    return SX(Sparsity::scalar(), *this, false);
  }

  SXElement SXElement::zz_min(const SXElement &b) const {
    return BinarySX::create(OP_FMIN, *this, b);
  }

  SXElement SXElement::zz_max(const SXElement &b) const {
    return BinarySX::create(OP_FMAX, *this, b);
  }

  SXElement SXElement::zz_atan2(const SXElement &b) const {
    return BinarySX::create(OP_ATAN2, *this, b);
  }

  SXElement SXElement::printme(const SXElement &b) const {
    return BinarySX::create(OP_PRINTME, *this, b);
  }

  SXElement SXElement::zz_power(const SXElement& n) const {
    if (n->isConstant()) {
      if (n->isInteger()) {
        int nn = n->getIntValue();
        if (nn == 0) {
          return 1;
        } else if (nn>100 || nn<-100) { // maximum depth
          return BinarySX::create(OP_CONSTPOW, *this, nn);
        } else if (nn<0) { // negative power
          return 1/pow(*this, -nn);
        } else if (nn%2 == 1) { // odd power
          return *this*pow(*this, nn-1);
        } else { // even power
          SXElement rt = pow(*this, nn/2);
          return rt*rt;
        }
      } else if (n->getValue()==0.5) {
        return sqrt(*this);
      } else {
        return BinarySX::create(OP_CONSTPOW, *this, n);
      }
    } else {
      return BinarySX::create(OP_POW, *this, n);
    }
  }

  SXElement SXElement::__constpow__(const SXElement& n) const {
    return BinarySX::create(OP_CONSTPOW, *this, n);
  }

  SXElement SXElement::constpow(const SXElement& n) const {
    return BinarySX::create(OP_CONSTPOW, *this, n);
  }

  SXElement SXElement::zz_not() const {
    if (isOp(OP_NOT)) {
      return getDep();
    } else {
      return UnarySX::create(OP_NOT, *this);
    }
  }

  SXElement SXElement::zz_and(const SXElement& y) const {
    return BinarySX::create(OP_AND, *this, y);
  }

  SXElement SXElement::zz_or(const SXElement& y) const {
    return BinarySX::create(OP_OR, *this, y);
  }

  SXElement SXElement::zz_if_else_zero(const SXElement& y) const {
    if (y->isZero()) {
      return y;
    } else if (isConstant()) {
      if (getValue()!=0) return y;
      else              return 0;
    } else {
      return BinarySX::create(OP_IF_ELSE_ZERO, *this, y);
    }
  }

  int SXElement::getTemp() const {
    return (*this)->temp;
  }

  void SXElement::setTemp(int t) {
    (*this)->temp = t;
  }

  bool SXElement::marked() const {
    return (*this)->marked();
  }

  void SXElement::mark() {
    (*this)->mark();
  }

  template<>
  void SX::setMaxNumCallsInPrint(long num) {
    SXNode::max_num_calls_in_print_ = num;
  }

  template<>
  long SX::getMaxNumCallsInPrint() {
    return SXNode::max_num_calls_in_print_;
  }

  template<>
  void SX::setEqualityCheckingDepth(int eq_depth) {
    SXNode::eq_depth_ = eq_depth;
  }

  template<>
  int SX::getEqualityCheckingDepth() {
    return SXNode::eq_depth_;
  }

  template<>
  SX GenericMatrix<SX>::sym(const std::string& name, const Sparsity& sp) {
    // Create a dense n-by-m matrix
    std::vector<SXElement> retv;

    // Check if individial names have been provided
    if (name[0]=='[') {

      // Make a copy of the string and modify it as to remove the special characters
      string modname = name;
      for (string::iterator it=modname.begin(); it!=modname.end(); ++it) {
        switch (*it) {
        case '(': case ')': case '[': case ']': case '{': case '}': case ',': case ';': *it = ' ';
        }
      }

      istringstream iss(modname);
      string varname;

      // Loop over elements
      while (!iss.fail()) {
        // Read the name
        iss >> varname;

        // Append to the return vector
        if (!iss.fail())
          retv.push_back(SXElement::sym(varname));
      }
    } else if (sp.isScalar(true)) {
      retv.push_back(SXElement::sym(name));
    } else {
      // Scalar
      std::stringstream ss;
      for (int k=0; k<sp.nnz(); ++k) {
        ss.str("");
        ss << name << "_" << k;
        retv.push_back(SXElement::sym(ss.str()));
      }
    }

    // Determine dimensions automatically if empty
    if (sp.isScalar(true)) {
      return SX(retv);
    } else {
      return SX(sp, retv, false);
    }
  }

  bool SXElement::isRegular() const {
    if (isConstant()) {
      return !(isNan() || isInf() || isMinusInf());
    } else {
      casadi_error("Cannot check regularity for symbolic SXElement");
    }
  }

  template<>
  bool SX::isRegular() const {
    // First pass: ignore symbolics
    for (int i=0; i<nnz(); ++i) {
      const SXElement& x = at(i);
      if (x.isConstant()) {
        if (x.isNan() || x.isInf() || x.isMinusInf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (int i=0; i<nnz(); ++i) {
      if (!at(i).isRegular()) return false;
    }
    return true;
  }

  template<>
  bool SX::isSmooth() const {
    // Make a function
    SXFunction temp(SX(), *this);
    temp.init();

    // Run the function on the temporary variable
    return temp->isSmooth();
  }

  template<>
  long SX::getElementHash() const {
    return toScalar().__hash__();
  }

  template<>
  bool SX::isLeaf() const {
    return toScalar().isLeaf();
  }

  template<>
  bool SX::isCommutative() const {
    return toScalar().isCommutative();
  }

  template<>
  bool SX::isSymbolic() const {
    if (isDense()) {
      return isValidInput();
    } else {
      return false;
    }
  }

  template<>
  bool SX::isValidInput() const {
    for (int k=0; k<nnz(); ++k) // loop over non-zero elements
      if (!at(k)->isSymbolic()) // if an element is not symbolic
        return false;

    return true;
  }

  template<> bool SX::hasDuplicates() {
    bool has_duplicates = false;
    for (vector<SXElement>::iterator it = begin(); it != end(); ++it) {
      bool is_duplicate = it->getTemp()!=0;
      if (is_duplicate) {
        cerr << "Duplicate expression: " << *it << endl;
      }
      has_duplicates = has_duplicates || is_duplicate;
      it->setTemp(1);
    }
    return has_duplicates;
  }

  template<> void SX::resetInput() {
    for (vector<SXElement>::iterator it = begin(); it != end(); ++it) {
      it->setTemp(0);
    }
  }

  template<>
  double SX::getValue(int k) const {
    return at(k).getValue();
  }

  template<>
  int SX::getIntValue() const {
    return toScalar().getIntValue();
  }

  template<>
  std::vector<double> SX::nonzeros() const {
    std::vector<double> ret(nnz());
    for (size_t i=0; i<ret.size(); ++i) {
      ret[i] = at(i).getValue();
    }
    return ret;
  }

  template<>
  std::vector<int> SX::nonzeros_int() const {
    std::vector<int> ret(nnz());
    for (size_t i=0; i<ret.size(); ++i) {
      ret[i] = at(i).getIntValue();
    }
    return ret;
  }

  template<>
  std::string SX::getName() const {
    return toScalar().getName();
  }

  template<>
  SX SX::getDep(int ch) const {
    return toScalar().getDep(ch);
  }

  template<>
  int SX::getNdeps() const {
    return toScalar().getNdeps();
  }

  template<>
  void SX::zz_expand(SX &ww, SX& tt) const {
    const SX& ex2 = *this;
    casadi_assert(ex2.isScalar());
    SXElement ex = ex2.toScalar();

    // Terms, weights and indices of the nodes that are already expanded
    std::vector<std::vector<SXNode*> > terms;
    std::vector<std::vector<double> > weights;
    std::map<SXNode*, int> indices;

    // Stack of nodes that are not yet expanded
    std::stack<SXNode*> to_be_expanded;
    to_be_expanded.push(ex.get());

    while (!to_be_expanded.empty()) { // as long as there are nodes to be expanded

      // Check if the last element on the stack is already expanded
      if (indices.find(to_be_expanded.top()) != indices.end()) {
        // Remove from stack
        to_be_expanded.pop();
        continue;
      }

      // Weights and terms
      std::vector<double> w; // weights
      std::vector<SXNode*> f; // terms

      if (to_be_expanded.top()->isConstant()) { // constant nodes are seen as multiples of one
        w.push_back(to_be_expanded.top()->getValue());
        f.push_back(casadi_limits<SXElement>::one.get());
      } else if (to_be_expanded.top()->isSymbolic()) {
        // symbolic nodes have weight one and itself as factor
        w.push_back(1);
        f.push_back(to_be_expanded.top());
      } else { // binary node

        casadi_assert(to_be_expanded.top()->hasDep()); // make sure that the node is binary

        // Check if addition, subtracton or multiplication
        SXNode* node = to_be_expanded.top();
        // If we have a binary node that we can factorize
        if (node->getOp() == OP_ADD || node->getOp() == OP_SUB ||
           (node->getOp() == OP_MUL  && (node->dep(0)->isConstant() ||
                                         node->dep(1)->isConstant()))) {
          // Make sure that both children are factorized, if not - add to stack
          if (indices.find(node->dep(0).get()) == indices.end()) {
            to_be_expanded.push(node->dep(0).get());
            continue;
          }
          if (indices.find(node->dep(1).get()) == indices.end()) {
            to_be_expanded.push(node->dep(1).get());
            continue;
          }

          // Get indices of children
          int ind1 = indices[node->dep(0).get()];
          int ind2 = indices[node->dep(1).get()];

          // If multiplication
          if (node->getOp() == OP_MUL) {
            double fac;
            if (node->dep(0)->isConstant()) { // Multiplication where the first factor is a constant
              fac = node->dep(0)->getValue();
              f = terms[ind2];
              w = weights[ind2];
            } else { // Multiplication where the second factor is a constant
              fac = node->dep(1)->getValue();
              f = terms[ind1];
              w = weights[ind1];
            }
            for (int i=0; i<w.size(); ++i) w[i] *= fac;

          } else { // if addition or subtraction
            if (node->getOp() == OP_ADD) {          // Addition: join both sums
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];    w.insert(w.end(), weights[ind2].begin(), weights[ind2].end());
            } else {      // Subtraction: join both sums with negative weights for second term
              f = terms[ind1];      f.insert(f.end(), terms[ind2].begin(), terms[ind2].end());
              w = weights[ind1];
              w.reserve(f.size());
              for (int i=0; i<weights[ind2].size(); ++i) w.push_back(-weights[ind2][i]);
            }
            // Eliminate multiple elements
            std::vector<double> w_new; w_new.reserve(w.size());   // weights
            std::vector<SXNode*> f_new;  f_new.reserve(f.size());   // terms
            std::map<SXNode*, int> f_ind; // index in f_new

            for (int i=0; i<w.size(); i++) {
              // Try to locate the node
              std::map<SXNode*, int>::iterator it = f_ind.find(f[i]);
              if (it == f_ind.end()) { // if the term wasn't found
                w_new.push_back(w[i]);
                f_new.push_back(f[i]);
                f_ind[f[i]] = f_new.size()-1;
              } else { // if the term already exists
                w_new[it->second] += w[i]; // just add the weight
              }
            }
            w = w_new;
            f = f_new;
          }
        } else { // if we have a binary node that we cannot factorize
          // By default,
          w.push_back(1);
          f.push_back(node);

        }
      }

      // Save factorization of the node
      weights.push_back(w);
      terms.push_back(f);
      indices[to_be_expanded.top()] = terms.size()-1;

      // Remove node from stack
      to_be_expanded.pop();
    }

    // Save expansion to output
    int thisind = indices[ex.get()];
    ww = SX(weights[thisind]);

    vector<SXElement> termsv(terms[thisind].size());
    for (int i=0; i<termsv.size(); ++i)
      termsv[i] = SXElement::create(terms[thisind][i]);
    tt = SX(termsv);
  }

  template<>
  SX SX::zz_pw_const(const SX &tval,
                     const SX &val) const {
    const SX &t = *this;
    // number of intervals
    int n = val.numel();

    casadi_assert_message(t.isScalar(), "t must be a scalar");
    casadi_assert_message(tval.numel() == n-1, "dimensions do not match");

    SX ret = val.at(0);
    for (int i=0; i<n-1; ++i) {
      ret += (val(0, i+1)-val(0, i)) * (t>=tval(0, i));
    }

    return ret;
  }

  template<>
  SX SX::zz_pw_lin(const SX &tval,
                   const SX &val) const {
    const SX &t = *this;

    // Number of points
    int N = tval.numel();
    casadi_assert_message(N>=2, "pw_lin: N>=2");

    // Gradient for each line segment
    SX g = SX(1, N-1);
    for (int i=0; i<N-1; ++i) {
      g(0, i) = (val(0, i+1)- val(0, i))/(tval(0, i+1)-tval(0, i));
    }

    // Line segments
    SX lseg = SX(1, N-1);
    for (int i=0; i<N-1; ++i)
      lseg(0, i) = val(0, i) + g(0, i)*(t-tval(0, i));

    // interior time points
    SX tint = tval(0, range(N-2));

    // Return piecewise linear function
    return pw_const(t, tint, lseg);
  }

  template<>
  SX SX::zz_if_else(const SX &if_true,
                    const SX &if_false) const {
    return if_else_zero(*this, if_true) + if_else_zero(!*this, if_false);
  }

  template<>
  SX SX::zz_heaviside() const {
    return (1+sign(*this))/2;
  }

  template<>
  SX SX::zz_rectangle() const {
    return 0.5*(sign(*this+0.5)-sign(*this-0.5));
  }

  template<>
  SX SX::zz_triangle() const {
    return rectangle(toScalar()/2)*(1-abs(toScalar()));
  }

  template<>
  SX SX::zz_ramp() const {
    return *this*heaviside(*this);
  }

  template<>
  SX SX::zz_gauss_quadrature(const SX &x,
                             const SX &a,
                             const SX &b, int order,
                             const SX& w) const {
    const SX &f = *this;
    casadi_assert_message(order == 5, "gauss_quadrature: order must be 5");
    casadi_assert_message(w.isEmpty(), "gauss_quadrature: empty weights");

    // Change variables to [-1, 1]
    if (!isEqual(a.toScalar(), -1) || !isEqual(b.toScalar(), 1)) {
      SX q1 = (b-a)/2;
      SX q2 = (b+a)/2;

      SXFunction fcn(x, f);
      fcn.init();

      return q1*gauss_quadrature(fcn(q1*x+q2).at(0), x, -1, 1);
    }

    // Gauss points
    vector<double> xi;
    xi.push_back(-sqrt(5 + 2*sqrt(10.0/7))/3);
    xi.push_back(-sqrt(5 - 2*sqrt(10.0/7))/3);
    xi.push_back(0);
    xi.push_back(sqrt(5 - 2*sqrt(10.0/7))/3);
    xi.push_back(sqrt(5 + 2*sqrt(10.0/7))/3);

    // Gauss weights
    vector<double> wi;
    wi.push_back((322-13*sqrt(70.0))/900.0);
    wi.push_back((322+13*sqrt(70.0))/900.0);
    wi.push_back(128/225.0);
    wi.push_back((322+13*sqrt(70.0))/900.0);
    wi.push_back((322-13*sqrt(70.0))/900.0);

    // Evaluate at the Gauss points
    SXFunction fcn(x, f);
    vector<SXElement> f_val(5);
    for (int i=0; i<5; ++i)
      f_val[i] = fcn(SX(xi[i])).at(0).toScalar();

    // Weighted sum
    SXElement sum;
    for (int i=0; i<5; ++i)
      sum += wi[i]*f_val[i];

    return sum;
  }

  template<>
  SX SX::zz_simplify() const {
    SX ex = *this;
    for (int el=0; el<ex.nnz(); ++el) ex.at(el) = simplify(ex.at(el));
    return ex;
  }

  template<>
  SX SX::zz_substitute(const SX& v,
                       const SX& vdef) const {
    return substitute(vector<SX>(1, *this), vector<SX>(1, v), vector<SX>(1, vdef)).front();
  }

  template<>
  std::vector<SX >
  SX::zz_substitute(const std::vector<SX >& ex,
                    const std::vector<SX >& v,
                    const std::vector<SX >& vdef) {

    // Assert consistent dimensions
    casadi_assert_warning(v.size()==vdef.size(), "subtitute: number of symbols to replace ( "
                          << v.size() << ") must match number of expressions (" << vdef.size()
                          << ") to replace them with.");

    // Quick return if all equal
    bool all_equal = true;
    for (int k=0; k<v.size(); ++k) {
      if (v[k].shape()!=vdef[k].shape() || !isEqual(v[k], vdef[k])) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) return ex;

    // Check sparsities
    for (int k=0; k<v.size(); ++k) {
      if (v[k].sparsity()!=vdef[k].sparsity()) {
        // Expand vdef to sparsity of v if vdef is scalar
        if (vdef[k].isScalar() && vdef[k].nnz()==1) {
          std::vector<SX> vdef_mod = vdef;
          vdef_mod[k] = SX(v[k].sparsity(), vdef[k].at(0), false);
          return substitute(ex, v, vdef_mod);
        } else {
          casadi_error("subsitute(ex, v, vdef): sparsities of v and vdef must match. Got v: "
                       << v[k].dimString() << " and " << "vdef: " << vdef[k].dimString() << ".");
        }
      }
    }


    // Otherwise, evaluate symbolically
    SXFunction F(v, ex);
    F.init();
    return F(vdef);
  }

  template<>
  void SX::zz_substituteInPlace(const std::vector<SX >& v,
                                std::vector<SX >& vdef,
                                std::vector<SX >& ex,
                                bool reverse) {
    // Assert correctness
    casadi_assert(v.size()==vdef.size());
    for (int i=0; i<v.size(); ++i) {
      casadi_assert_message(v[i].isSymbolic(), "the variable is not symbolic");
      casadi_assert_message(v[i].sparsity() == vdef[i].sparsity(), "the sparsity patterns of the "
                            "expression and its defining bexpression do not match");
    }

    // Quick return if empty or single expression
    if (v.empty()) return;

    // Function inputs
    std::vector<SX> f_in;
    if (!reverse) f_in.insert(f_in.end(), v.begin(), v.end());

    // Function outputs
    std::vector<SX> f_out = vdef;
    f_out.insert(f_out.end(), ex.begin(), ex.end());

    // Write the mapping function
    SXFunction f(f_in, f_out);
    f.init();

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = f.algorithm();
    vector<SXElement> work(f.getWorkSize());

    // Iterator to the binary operations
    vector<SXElement>::const_iterator b_it=f->operations_.begin();

    // Iterator to stack of constants
    vector<SXElement>::const_iterator c_it = f->constants_.begin();

    // Iterator to free variables
    vector<SXElement>::const_iterator p_it = f->free_vars_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_INPUT:
        // reverse is false, substitute out
        work[it->i0] = vdef.at(it->i1).at(it->i2);
        break;
      case OP_OUTPUT:
        if (it->i0 < v.size()) {
          vdef.at(it->i0).at(it->i2) = work[it->i1];
          if (reverse) {
            // Use the new variable henceforth, substitute in
            work[it->i1] = v.at(it->i0).at(it->i2);
          }
        } else {
          // Auxillary output
          ex.at(it->i0 - v.size()).at(it->i2) = work[it->i1];
        }
        break;
      case OP_CONST:      work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }

          // Avoid creating duplicates
          const int depth = 2; // NOTE: a higher depth could possibly give more savings
          work[it->i0].assignIfDuplicate(*b_it++, depth);
        }
      }
    }
  }

  template<>
  SX SX::zz_spy() const {
    SX s = SX(size1(), size2());
    for (int i=0; i<size2(); ++i)
      for (int j=0; j<size1(); ++j)
        if (!(*this)(j, i).toScalar()->isZero())
          s(j, i) = 1;
    return s;
  }

  template<>
  bool SX::zz_dependsOn(const SX &arg) const {
    const SX& ex = *this;
    if (ex.nnz()==0) return false;

    // Construct a temporary algorithm
    SXFunction temp(arg, ex);
    temp.init();
    temp.spInit(true);

    bvec_t* input_ =  get_bvec_t(temp.input().data());
    // Make a column with all variables active
    std::fill(input_, input_+temp.input().nnz(), bvec_t(1));
    bvec_t* output_ = get_bvec_t(temp.output().data());
    // Perform a single dependency sweep
    temp.spEvaluate(true);

    // Loop over results
    for (int i=0;i<temp.output().nnz();++i) {
      if (output_[i]) return true;
    }

    return false;
  }

  template<>
  SX SX::zz_jacobian(const SX &arg) const {
    SXFunction temp(arg, *this); // make a runtime
    temp.init();
    return temp.jac();
  }

  template<>
  SX SX::zz_gradient(const SX &arg) const {
    SXFunction temp(arg, *this); // make a runtime
    temp.init();
    return temp.grad();
  }

  template<>
  SX SX::zz_tangent(const SX &arg) const {
    SXFunction temp(arg, *this); // make a runtime
    temp.init();
    return temp.tang();
  }

  template<>
  SX SX::zz_hessian(const SX &arg) const {
    SX H, g;
    hessian(*this, arg, H, g);
    return H;
  }

  template<>
  void SX::zz_hessian(const SX &arg, SX &H, SX &g) const {
    g = gradient(*this, arg);

    SXFunction temp(arg, g); // make a runtime
    temp.init();
    H = temp.jac(0, 0, false, true);
  }

  template<>
  SX SX::zz_jacobianTimesVector(const SX &arg, const SX &v, bool transpose_jacobian) const {
    SXFunction f(arg, *this);
    f.init();

    // Split up v
    vector<SX> vv = horzsplit(v);

    // Make sure well-posed
    casadi_assert(vv.size() >= 1);
    casadi_assert(isVector());
    casadi_assert(arg.isVector());
    if (transpose_jacobian) {
      casadi_assert(v.size1()==size1());
    } else {
      casadi_assert(v.size1()==arg.size1());
    }

    // Number of sensitivities
    int nfsens = transpose_jacobian ? 0 : vv.size();
    int nasens = transpose_jacobian ? vv.size() : 0;

    // Assemble arguments and directional derivatives
    vector<SX> argv = f.inputExpr();
    vector<SX> resv = f.outputExpr();
    vector<vector<SX> > fseed(nfsens, argv), fsens(nfsens, resv),
        aseed(nasens, resv), asens(nasens, argv);

    for (int dir=0; dir<vv.size(); ++dir) {
      if (transpose_jacobian) {
        aseed[dir][0].set(vv[dir]);
      } else {
        fseed[dir][0].set(vv[dir]);
      }
    }

    // Evaluate with directional derivatives, output is the same as the funciton inputs
    f.callDerivative(argv, resv, fseed, fsens, aseed, asens);

    // Get the results
    for (int dir=0; dir<vv.size(); ++dir) {
      if (transpose_jacobian) {
        vv[dir] = asens[dir][0];
      } else {
        vv[dir] = fsens[dir][0];
      }
    }
    return horzcat(vv);
  }

  template<>
  SX SX::zz_taylor(const SX& x,
                   const SX& a, int order) const {
    casadi_assert(x.isScalar() && a.isScalar());
    if (nnz()!=numel())
      throw CasadiException("taylor: not implemented for sparse matrices");
    SX ff = vec(T());

    SX result = substitute(ff, x, a);
    double nf=1;
    SX dx = (x-a);
    SX dxa = (x-a);
    for (int i=1;i<=order;i++) {
      ff = jacobian(ff, x);
      nf*=i;
      result+=1/nf * substitute(ff, x, a) * dxa;
      dxa*=dx;
    }
    return reshape(result, size2(), size1()).T();
  }

  template<>
  SX SX::zz_mtaylor(const SX& x, const SX& a, int order) const {
    return mtaylor(*this, x, a, order, std::vector<int>(x.nnz(), 1));
  }

  SX mtaylor_recursive(const SX& ex, const SX& x, const SX& a, int order,
                       const std::vector<int>&order_contributions,
                       const SXElement & current_dx=casadi_limits<SXElement>::one,
                       double current_denom=1, int current_order=1) {
    SX result = substitute(ex, x, a)*current_dx/current_denom;
    for (int i=0;i<x.nnz();i++) {
      if (order_contributions[i]<=order) {
        result += mtaylor_recursive(
                                    jacobian(ex, x.at(i)),
                                    x, a,
                                    order-order_contributions[i],
                                    order_contributions,
                                    current_dx*(x.at(i)-a.at(i)),
                                    current_denom*current_order, current_order+1);
      }
    }
    return result;
  }

  template<>
  SX SX::zz_mtaylor(const SX& x, const SX& a, int order,
                    const std::vector<int>& order_contributions) const {
    casadi_assert_message(nnz()==numel() && x.nnz()==x.numel(),
                          "mtaylor: not implemented for sparse matrices");

    casadi_assert_message(x.nnz()==order_contributions.size(),
                          "mtaylor: number of non-zero elements in x (" <<  x.nnz()
                          << ") must match size of order_contributions ("
                          << order_contributions.size() << ")");

    return reshape(mtaylor_recursive(vec(*this), x, a, order,
                                     order_contributions),
                   size2(), size1()).T();
  }

  template<>
  int SX::zz_countNodes() const {
    SXFunction f(SX(), *this);
    f.init();
    return f.countNodes();
  }

  template<>
  std::string
  SX::zz_getOperatorRepresentation(const std::vector<std::string>& args) const {
    SXElement x = toScalar();
    if (!x.hasDep())
        throw CasadiException("getOperatorRepresentation: SXElement must be binary operator");
    if (args.size() == 0 || (casadi_math<double>::ndeps(x.getOp())==2 && args.size() < 2))
        throw CasadiException("getOperatorRepresentation: not enough arguments supplied");
    std::stringstream s;
    casadi_math<double>::print(x.getOp(), s, args[0], args[1]);
    return s.str();
  }

  template<>
  std::vector<SX> SX::zz_getSymbols() const {
    SXFunction f(std::vector<SX>(), *this);
    f.init();
    return std::vector<SX>(1, f.getFree());
  }

  template<>
  std::vector<SX >
  SX::zz_getSymbols(const std::vector<SX >& e) {
    throw CasadiException("\"getSymbols\" not defined for instantiation");
    return std::vector<SX >();
  }

  template<>
  void SX::zz_extractShared(std::vector<SX >& ex,
                            std::vector<SX >& v_sx,
                            std::vector<SX >& vdef_sx,
                            const std::string& v_prefix,
                            const std::string& v_suffix) {

    // Sort the expression
    SXFunction f(vector<SX>(), ex);
    f.init();

    // Get references to the internal data structures
    const vector<ScalarAtomic>& algorithm = f.algorithm();
    vector<SXElement> work(f.getWorkSize());
    vector<SXElement> work2 = work;

    // Iterator to the binary operations
    vector<SXElement>::const_iterator b_it=f->operations_.begin();

    // Iterator to stack of constants
    vector<SXElement>::const_iterator c_it = f->constants_.begin();

    // Iterator to free variables
    vector<SXElement>::const_iterator p_it = f->free_vars_.begin();

    // Count how many times an expression has been used
    vector<int> usecount(work.size(), 0);

    // Evaluate the algorithm
    vector<SXElement> v, vdef;
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      // Increase usage counters
      switch (it->op) {
      case OP_CONST:
      case OP_PARAMETER:
        break;
        CASADI_MATH_BINARY_BUILTIN // Binary operation
          if (usecount[it->i2]==0) {
            usecount[it->i2]=1;
          } else if (usecount[it->i2]==1) {
            // Get a suitable name
            vdef.push_back(work[it->i2]);
            usecount[it->i2]=-1; // Extracted, do not extract again
          }
        // fall-through
      case OP_OUTPUT:
      default: // Unary operation, binary operation or output
        if (usecount[it->i1]==0) {
          usecount[it->i1]=1;
        } else if (usecount[it->i1]==1) {
          vdef.push_back(work[it->i1]);
          usecount[it->i1]=-1; // Extracted, do not extract again
        }
      }

      // Perform the operation
      switch (it->op) {
      case OP_OUTPUT:
        break;
      case OP_CONST:
      case OP_PARAMETER:
        usecount[it->i0] = -1; // Never extract since it is a primitive type
        break;
      default:
        work[it->i0] = *b_it++;
        usecount[it->i0] = 0; // Not (yet) extracted
        break;
      }
    }

    // Create intermediate variables
    stringstream v_name;
    for (int i=0; i<vdef.size(); ++i) {
      v_name.str(string());
      v_name << v_prefix << i << v_suffix;
      v.push_back(SXElement::sym(v_name.str()));
    }

    // Mark the above expressions
    for (int i=0; i<vdef.size(); ++i) {
      vdef[i].setTemp(i+1);
    }

    // Save the marked nodes for later cleanup
    vector<SXElement> marked = vdef;

    // Reset iterator
    b_it=f->operations_.begin();

    // Evaluate the algorithm
    for (vector<ScalarAtomic>::const_iterator it=algorithm.begin(); it<algorithm.end(); ++it) {
      switch (it->op) {
      case OP_OUTPUT:     ex.at(it->i0).at(it->i2) = work[it->i1];      break;
      case OP_CONST:      work2[it->i0] = work[it->i0] = *c_it++; break;
      case OP_PARAMETER:  work2[it->i0] = work[it->i0] = *p_it++; break;
      default:
        {
          switch (it->op) {
            CASADI_MATH_FUN_BUILTIN(work[it->i1], work[it->i2], work[it->i0])
              }
          work2[it->i0] = *b_it++;

          // Replace with intermediate variables
          int ind = work2[it->i0].getTemp()-1;
          if (ind>=0) {
            vdef.at(ind) = work[it->i0];
            work[it->i0] = v.at(ind);
          }
        }
      }
    }

    // Unmark the expressions
    for (vector<SXElement>::iterator it=marked.begin(); it!=marked.end(); ++it) {
      it->setTemp(0);
    }

    // Save v, vdef
    v_sx.resize(v.size());
    std::copy(v.begin(), v.end(), v_sx.begin());
    vdef_sx.resize(vdef.size());
    std::copy(vdef.begin(), vdef.end(), vdef_sx.begin());
  }

  template<>
  void SX::zz_printCompact(std::ostream &stream) const {
    // Extract shared subexpressions from ex
    vector<SX> v, vdef;
    vector<SX> ex_extracted(1, *this);
    extractShared(ex_extracted, v, vdef, "@", "");

    // Print the expression without shared subexpressions
    ex_extracted.at(0).print(stream);

    // Print the shared subexpressions
    if (!v.empty()) {
      stream << endl << "where:" << endl;
      for (int i=0; i<v.size(); ++i) {
        stream << v[i].toScalar() << " := " << vdef[i].toScalar() << endl;
      }
    }
  }

  template<>
  SX SX::zz_poly_coeff(const SX&x) const {
    casadi_assert(isScalar());
    casadi_assert(x.isScalar());
    casadi_assert(x.isSymbolic());

    SX ret;

    SXFunction f(x, *this);
    f.init();
    int mult = 1;
    bool success = false;
    for (int i=0;i<1000;++i) {
      ret.append(f(SX(0)).at(0)/mult);
      SX j = f.jac();
      if (j.nnz()==0) {
        success = true;
        break;
      }
      f = SXFunction(x, j);
      f.init();
      mult*=i+1;
    }

    if (!success) casadi_error("poly: supplied expression does not appear to be polynomial.");

    std::reverse(ret.data().begin(), ret.data().end());

    return ret;
  }

  template<>
  SX SX::zz_poly_roots() const {
    const SX& p = *this;
    casadi_assert_message(p.size2()==1,
                          "poly_root(): supplied parameter must be column vector but got "
                          << p.dimString() << ".");
    casadi_assert(p.isDense());
    if (p.size1()==2) { // a*x + b
      SX a = p(0);
      SX b = p(1);
      return -b/a;
    } else if (p.size1()==3) { // a*x^2 + b*x + c
      SX a = p(0);
      SX b = p(1);
      SX c = p(2);
      SX ds = sqrt(b*b-4*a*c);
      SX bm = -b;
      SX a2 = 2*a;
      SX ret;
      ret.append((bm-ds)/a2);
      ret.append((bm+ds)/a2);
      return ret;
    } else if (p.size1()==4) {
      // www.cs.iastate.edu/~cs577/handouts/polyroots.pdf
      SX ai = 1/p(0);

      SX p_ = p(1)*ai;
      SX q  = p(2)*ai;
      SX r  = p(3)*ai;

      SX pp = p_*p_;

      SX a = q - pp/3;
      SX b = r + 2.0/27*pp*p_-p_*q/3;

      SX a3 = a/3;

      SX phi = acos(-b/2/sqrt(-a3*a3*a3));

      SX ret;
      ret.append(cos(phi/3));
      ret.append(cos((phi+2*M_PI)/3));
      ret.append(cos((phi+4*M_PI)/3));
      ret*= 2*sqrt(-a3);

      ret-= p_/3;
      return ret;
    } else if (p.size1()==5) {
      SX ai = 1/p(0);
      SX b = p(1)*ai;
      SX c = p(2)*ai;
      SX d = p(3)*ai;
      SX e = p(4)*ai;

      SX bb= b*b;
      SX f = c - (3*bb/8);
      SX g = d + (bb*b / 8) - b*c/2;
      SX h = e - (3*bb*bb/256) + (bb * c/16) - (b*d/4);
      SX poly;
      poly.append(1);
      poly.append(f/2);
      poly.append((f*f -4*h)/16);
      poly.append(-g*g/64);
      SX y = poly_roots(poly);

      SX r0 = y(0);
      SX r1 = y(2);

      SX p = sqrt(r0); // two non-zero-roots
      SX q = sqrt(r1);

      SX r = -g/(8*p*q);

      SX s = b/4;

      SX ret;
      ret.append(p + q + r -s);
      ret.append(p - q - r -s);
      ret.append(-p + q - r -s);
      ret.append(-p - q + r -s);

      return ret;
    } else if (isEqual(p(p.nnz()-1).at(0), 0)) {
      SX ret = poly_roots(p(range(p.nnz()-1)));
      ret.append(0);
      return ret;
    } else {
      casadi_error("poly_root(): can only solve cases for first or second order polynomial. "
                   "Got order " << p.size1()-1 << ".");
    }

  }

  template<>
  SX SX::zz_eig_symbolic() const {
    const SX& m = *this;
    casadi_assert_message(m.size1()==m.size2(), "eig(): supplied matrix must be square");

    SX ret;

    /// Bring m in block diagonal form, calculating eigenvalues of each block separately
    std::vector<int> offset;
    std::vector<int> index;
    int nb = m.sparsity().stronglyConnectedComponents(offset, index);

    SX m_perm = m(offset, offset);

    SX l = SX::sym("l");

    for (int k=0;k<nb;++k) {
      std::vector<int> r = range(index.at(k), index.at(k+1));
      // det(lambda*I-m) = 0
      ret.append(poly_roots(poly_coeff(det(SX::eye(r.size())*l-m_perm(r, r)), l)));
    }

    return ret;
  }

  SXElement SXElement::zz_simplify() const {
    // Start by expanding the node to a weighted sum
    SX terms, weights;
    expand(*this, weights, terms);

    // Make a scalar product to get the simplified expression
    SX s = mul(terms.T(), weights);
    return s.toScalar();
  }

} // namespace casadi

using namespace casadi;
namespace std {
/**
  const bool numeric_limits<casadi::SXElement>::is_specialized = true;
  const int  numeric_limits<casadi::SXElement>::digits = 0;
  const int  numeric_limits<casadi::SXElement>::digits10 = 0;
  const bool numeric_limits<casadi::SXElement>::is_signed = false;
  const bool numeric_limits<casadi::SXElement>::is_integer = false;
  const bool numeric_limits<casadi::SXElement>::is_exact = false;
  const int numeric_limits<casadi::SXElement>::radix = 0;
  const int  numeric_limits<casadi::SXElement>::min_exponent = 0;
  const int  numeric_limits<casadi::SXElement>::min_exponent10 = 0;
  const int  numeric_limits<casadi::SXElement>::max_exponent = 0;
  const int  numeric_limits<casadi::SXElement>::max_exponent10 = 0;

  const bool numeric_limits<casadi::SXElement>::has_infinity = true;
  const bool numeric_limits<casadi::SXElement>::has_quiet_NaN = true;
  const bool numeric_limits<casadi::SXElement>::has_signaling_NaN = false;
  const float_denorm_style has_denorm = denorm absent;
  const bool numeric_limits<casadi::SXElement>::has_denorm_loss = false;

  const bool numeric_limits<casadi::SXElement>::is_iec559 = false;
  const bool numeric_limits<casadi::SXElement>::is_bounded = false;
  const bool numeric_limits<casadi::SXElement>::is_modulo = false;

  const bool numeric_limits<casadi::SXElement>::traps = false;
  const bool numeric_limits<casadi::SXElement>::tinyness_before = false;
  const float_round_style numeric_limits<casadi::SXElement>::round_style = round_toward_zero;
*/
  SXElement numeric_limits<SXElement>::infinity() throw() {
    return casadi::casadi_limits<SXElement>::inf;
  }

  SXElement numeric_limits<SXElement>::quiet_NaN() throw() {
    return casadi::casadi_limits<SXElement>::nan;
  }

  SXElement numeric_limits<SXElement>::min() throw() {
    return SXElement(numeric_limits<double>::min());
  }

  SXElement numeric_limits<SXElement>::max() throw() {
    return SXElement(numeric_limits<double>::max());
  }

  SXElement numeric_limits<SXElement>::epsilon() throw() {
    return SXElement(numeric_limits<double>::epsilon());
  }

  SXElement numeric_limits<SXElement>::round_error() throw() {
    return SXElement(numeric_limits<double>::round_error());
  }

} // namespace std

