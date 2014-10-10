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
    if (!isEqual(scalar, 0) && isEqual(scalar, depth)) {
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
    if (node->hasDep() && node->getOp() == OP_NEG)
      return node->dep(0);
    else if (node->isZero())
      return 0;
    else if (node->isMinusOne())
      return 1;
    else if (node->isOne())
      return -1;
    else
      return UnarySX::create(OP_NEG, *this);
  }

  SXElement SXElement::sign() const {
    return UnarySX::create(OP_SIGN, *this);
  }

  SXElement SXElement::__copysign__(const SXElement &y) const {
    return BinarySX::create(OP_COPYSIGN, *this, y);
  }

  SXElement SXElement::erfinv() const {
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

  SXElement SXElement::__add__(const SXElement& y) const {
    // NOTE: Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_ADD, *this, y);

    if (node->isZero())
      return y;
    else if (y->isZero()) // term2 is zero
      return *this;
    else if (y.hasDep() && y.getOp()==OP_NEG) // x + (-y) -> x - y
      return __sub__(-y);
    else if (hasDep() && getOp()==OP_NEG) // (-x) + y -> y - x
      return y.__sub__(getDep());
    else if (hasDep() && getOp()==OP_MUL &&
            y.hasDep() && y.getOp()==OP_MUL &&
            getDep(0).isConstant() && getDep(0).getValue()==0.5 &&
            y.getDep(0).isConstant() && y.getDep(0).getValue()==0.5 &&
            y.getDep(1).isEqual(getDep(1), SXNode::eq_depth_)) // 0.5x+0.5x = x
      return getDep(1);
    else if (hasDep() && getOp()==OP_DIV &&
            y.hasDep() && y.getOp()==OP_DIV &&
            getDep(1).isConstant() && getDep(1).getValue()==2 &&
            y.getDep(1).isConstant() && y.getDep(1).getValue()==2 &&
            y.getDep(0).isEqual(getDep(0), SXNode::eq_depth_)) // x/2+x/2 = x
      return getDep(0);
    else if (hasDep() && getOp()==OP_SUB && getDep(1).isEqual(y, SXNode::eq_depth_))
      return getDep(0);
    else if (y.hasDep() && y.getOp()==OP_SUB && isEqual(y.getDep(1), SXNode::eq_depth_))
      return y.getDep(0);
    else // create a new branch
      return BinarySX::create(OP_ADD, *this, y);
  }

  SXElement SXElement::__sub__(const SXElement& y) const {
    // Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_SUB, *this, y);

    if (y->isZero()) // term2 is zero
      return *this;
    if (node->isZero()) // term1 is zero
      return -y;
    if (isEqual(y, SXNode::eq_depth_)) // the terms are equal
      return 0;
    else if (y.hasDep() && y.getOp()==OP_NEG) // x - (-y) -> x + y
      return __add__(-y);
    else if (hasDep() && getOp()==OP_ADD && getDep(1).isEqual(y, SXNode::eq_depth_))
      return getDep(0);
    else if (hasDep() && getOp()==OP_ADD && getDep(0).isEqual(y, SXNode::eq_depth_))
      return getDep(1);
    else if (y.hasDep() && y.getOp()==OP_ADD && isEqual(y.getDep(1), SXNode::eq_depth_))
      return -y.getDep(0);
    else if (y.hasDep() && y.getOp()==OP_ADD && isEqual(y.getDep(0), SXNode::eq_depth_))
      return -y.getDep(1);
    else // create a new branch
      return BinarySX::create(OP_SUB, *this, y);
  }

  SXElement SXElement::__mul__(const SXElement& y) const {

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_MUL, *this, y);

    // Only simplifications that do not result in extra nodes area allowed
    if (y.isEqual(*this, SXNode::eq_depth_))
      return sq();
    else if (!isConstant() && y.isConstant())
      return y.__mul__(*this);
    else if (node->isZero() || y->isZero()) // one of the terms is zero
      return 0;
    else if (node->isOne()) // term1 is one
      return y;
    else if (y->isOne()) // term2 is one
      return *this;
    else if (y->isMinusOne())
      return -(*this);
    else if (node->isMinusOne())
      return -y;
    else if (y.hasDep() && y.getOp()==OP_INV)
      return (*this)/y.inv();
    else if (hasDep() && getOp()==OP_INV)
      return y/inv();
    else if (isConstant() && y.hasDep() && y.getOp()==OP_MUL && y.getDep(0).isConstant() &&
            getValue()*y.getDep(0).getValue()==1) // 5*(0.2*x) = x
      return y.getDep(1);
    else if (isConstant() && y.hasDep() && y.getOp()==OP_DIV && y.getDep(1).isConstant() &&
            getValue()==y.getDep(1).getValue()) // 5*(x/5) = x
      return y.getDep(0);
    else if (hasDep() && getOp()==OP_DIV && getDep(1).isEqual(y, SXNode::eq_depth_)) // ((2/x)*x)
      return getDep(0);
    else if (y.hasDep() && y.getOp()==OP_DIV &&
            y.getDep(1).isEqual(*this, SXNode::eq_depth_)) // ((2/x)*x)
      return y.getDep(0);
    else     // create a new branch
      return BinarySX::create(OP_MUL, *this, y);
  }


  bool SXElement::isDoubled() const {
    return isOp(OP_ADD) && node->dep(0).isEqual(node->dep(1), SXNode::eq_depth_);
  }

  SXElement SXElement::__div__(const SXElement& y) const {
    // Only simplifications that do not result in extra nodes area allowed

    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_DIV, *this, y);

    if (y->isZero()) // term2 is zero
      return casadi_limits<SXElement>::nan;
    else if (node->isZero()) // term1 is zero
      return 0;
    else if (y->isOne()) // term2 is one
      return *this;
    else if (y->isMinusOne())
      return -(*this);
    else if (isEqual(y, SXNode::eq_depth_)) // terms are equal
      return 1;
    else if (isDoubled() && y.isEqual(2))
      return node->dep(0);
    else if (isOp(OP_MUL) && y.isEqual(node->dep(0), SXNode::eq_depth_))
      return node->dep(1);
    else if (isOp(OP_MUL) && y.isEqual(node->dep(1), SXNode::eq_depth_))
      return node->dep(0);
    else if (node->isOne())
      return y.inv();
    else if (y.hasDep() && y.getOp()==OP_INV)
      return (*this)*y.inv();
    else if (isDoubled() && y.isDoubled())
      return node->dep(0) / y->dep(0);
    else if (y.isConstant() && hasDep() && getOp()==OP_DIV && getDep(1).isConstant() &&
            y.getValue()*getDep(1).getValue()==1) // (x/5)/0.2
      return getDep(0);
    else if (y.hasDep() && y.getOp()==OP_MUL &&
            y.getDep(1).isEqual(*this, SXNode::eq_depth_)) // x/(2*x) = 1/2
      return BinarySX::create(OP_DIV, 1, y.getDep(0));
    else if (hasDep() && getOp()==OP_NEG &&
            getDep(0).isEqual(y, SXNode::eq_depth_))      // (-x)/x = -1
      return -1;
    else if (y.hasDep() && y.getOp()==OP_NEG &&
            y.getDep(0).isEqual(*this, SXNode::eq_depth_))      // x/(-x) = 1
      return -1;
    else if (y.hasDep() && y.getOp()==OP_NEG && hasDep() &&
            getOp()==OP_NEG && getDep(0).isEqual(y.getDep(0), SXNode::eq_depth_))  // (-x)/(-x) = 1
      return 1;
    else if (isOp(OP_DIV) && y.isEqual(node->dep(0), SXNode::eq_depth_))
      return node->dep(1).inv();
    else // create a new branch
      return BinarySX::create(OP_DIV, *this, y);
  }

  SXElement SXElement::inv() const {
    if (node->hasDep() && node->getOp()==OP_INV) {
      return node->dep(0);
    } else {
      return UnarySX::create(OP_INV, *this);
    }
  }

  SX SXElement::fmin(const SX& b) const {
    return SX(*this).fmin(b);
  }
  SX SXElement::fmax(const SX& b) const {
    return SX(*this).fmax(b);
  }
  SX SXElement::constpow(const SX& n) const {
    return SX(*this).__constpow__(n);
  }
  SX SXElement::__copysign__(const SX& n) const {
    return SX(*this).__copysign__(n);
  }

  SX SXElement::arctan2(const SX& b) const {
    return SX(*this).arctan2(b);
  }

  SXElement SXElement::__le__(const SXElement& y) const {
    if ((y-(*this)).isNonNegative())
      return 1;
    else
      return BinarySX::create(OP_LE, *this, y);
  }

  SXElement SXElement::__lt__(const SXElement& y) const {
    if (((*this)-y).isNonNegative())
      return 0;
    else
      return BinarySX::create(OP_LT, *this, y);
  }

  SXElement SXElement::__eq__(const SXElement& y) const {
    if (isEqual(y))
      return 1;
    else
      return BinarySX::create(OP_EQ, *this, y);
  }

  SXElement SXElement::__ne__(const SXElement& y) const {
    if (isEqual(y))
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
    return node->isConstant() || node->isSymbolic();
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

  bool SXElement::isEqual(const SXElement& ex, int depth) const {
    if (node==ex.get())
      return true;
    else if (depth>0)
      return node->isEqual(ex.get(), depth);
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

  SXElement SXElement::exp() const {
    return UnarySX::create(OP_EXP, *this);
  }

  SXElement SXElement::log() const {
    return UnarySX::create(OP_LOG, *this);
  }

  SXElement SXElement::log10() const {
    return log()*(1/std::log(10.));
  }

  SXElement SXElement::sqrt() const {
    if (isOp(OP_SQ))
      return node->dep(0).fabs();
    else
      return UnarySX::create(OP_SQRT, *this);
  }

  SXElement SXElement::sq() const {
    if (isOp(OP_SQRT))
      return node->dep(0);
    else
      return UnarySX::create(OP_SQ, *this);
  }

  SXElement SXElement::sin() const {
    return UnarySX::create(OP_SIN, *this);
  }

  SXElement SXElement::cos() const {
    return UnarySX::create(OP_COS, *this);
  }

  SXElement SXElement::tan() const {
    return UnarySX::create(OP_TAN, *this);
  }

  SXElement SXElement::arcsin() const {
    return UnarySX::create(OP_ASIN, *this);
  }

  SXElement SXElement::arccos() const {
    return UnarySX::create(OP_ACOS, *this);
  }

  SXElement SXElement::arctan() const {
    return UnarySX::create(OP_ATAN, *this);
  }

  SXElement SXElement::sinh() const {
    if (node->isZero())
      return 0;
    else
      return UnarySX::create(OP_SINH, *this);
  }

  SXElement SXElement::cosh() const {
    if (node->isZero())
      return 1;
    else
      return UnarySX::create(OP_COSH, *this);
  }

  SXElement SXElement::tanh() const {
    if (node->isZero())
      return 0;
    else
      return UnarySX::create(OP_TANH, *this);
  }

  SXElement SXElement::arctanh() const {
    if (node->isZero())
      return 0;
    else
      return UnarySX::create(OP_ATANH, *this);
  }

  SXElement SXElement::arccosh() const {
    if (node->isOne())
      return 0;
    else
      return UnarySX::create(OP_ACOSH, *this);
  }

  SXElement SXElement::arcsinh() const {
    if (node->isZero())
      return 0;
    else
      return UnarySX::create(OP_ASINH, *this);
  }

  SXElement SXElement::floor() const {
    return UnarySX::create(OP_FLOOR, *this);
  }

  SXElement SXElement::ceil() const {
    return UnarySX::create(OP_CEIL, *this);
  }

  SXElement SXElement::fmod(const SXElement &b) const {
    return BinarySX::create(OP_FMOD, *this, b);
  }

  SXElement SXElement::erf() const {
    return UnarySX::create(OP_ERF, *this);
  }

  SXElement SXElement::fabs() const {
    if (isOp(OP_FABS) || isOp(OP_SQ))
      return *this;
    else
      return UnarySX::create(OP_FABS, *this);
  }

  SXElement::operator SX() const {
    return SX(Sparsity::scalar(), *this);
  }

  SXElement SXElement::fmin(const SXElement &b) const {
    return BinarySX::create(OP_FMIN, *this, b);
  }

  SXElement SXElement::fmax(const SXElement &b) const {
    return BinarySX::create(OP_FMAX, *this, b);
  }

  SXElement SXElement::arctan2(const SXElement &b) const {
    return BinarySX::create(OP_ATAN2, *this, b);
  }

  SXElement SXElement::printme(const SXElement &b) const {
    return BinarySX::create(OP_PRINTME, *this, b);
  }

  SXElement SXElement::__pow__(const SXElement& n) const {
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
        return sqrt();
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

  SXElement SXElement::logic_not() const {
    if (hasDep() && getOp() == OP_NOT) {
      return getDep();
    } else {
      return UnarySX::create(OP_NOT, *this);
    }
  }

  SXElement SXElement::logic_and(const SXElement& y) const {
    return BinarySX::create(OP_AND, *this, y);
  }

  SXElement SXElement::logic_or(const SXElement& y) const {
    return BinarySX::create(OP_OR, *this, y);
  }

  SXElement SXElement::if_else_zero(const SXElement& y) const {
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
      for (int k=0; k<sp.size(); ++k) {
        ss.str("");
        ss << name << "_" << k;
        retv.push_back(SXElement::sym(ss.str()));
      }
    }

    // Determine dimensions automatically if empty
    if (sp.isScalar(true)) {
      return SX(retv);
    } else {
      return SX(sp, retv);
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
    for (int i=0; i<size(); ++i) {
      const SXElement& x = at(i);
      if (x.isConstant()) {
        if (x.isNan() || x.isInf() || x.isMinusInf()) return false;
      }
    }
    // Second pass: don't ignore symbolics
    for (int i=0; i<size(); ++i) {
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
  bool SX::isSymbolic() const {
    if (isDense()) {
      return isSymbolicSparse();
    } else {
      return false;
    }
  }

  template<>
  bool SX::isSymbolicSparse() const {
    for (int k=0; k<size(); ++k) // loop over non-zero elements
      if (!at(k)->isSymbolic()) // if an element is not symbolic
        return false;

    return true;
  }

  template<>
  double SX::getValue() const {
    return toScalar().getValue();
  }

  template<>
  std::string SX::getName() const {
    return toScalar().getName();
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

