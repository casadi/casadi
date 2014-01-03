/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "sx.hpp"
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

using namespace std;
namespace CasADi{

  // Allocate storage for the caching
  CACHING_MAP<int,IntegerSX*> IntegerSX::cached_constants_;
  CACHING_MAP<double,RealtypeSX*> RealtypeSX::cached_constants_;

  SX::SX(){
    node = casadi_limits<SX>::nan.node;
    node->count++;
  }

  SX::SX(SXNode* node_, bool dummy) : node(node_){
    node->count++;
  }

  SX SX::create(SXNode* node){
    return SX(node,false);
  }

  SX::SX(const SX& scalar){
    node = scalar.node;
    node->count++;
  }

  SX::SX(double val){
    int intval = int(val);
    if(val-intval == 0){ // check if integer
      if(intval == 0)             node = casadi_limits<SX>::zero.node;
      else if(intval == 1)        node = casadi_limits<SX>::one.node;
      else if(intval == 2)        node = casadi_limits<SX>::two.node;
      else if(intval == -1)       node = casadi_limits<SX>::minus_one.node;
      else                        node = IntegerSX::create(intval);
      node->count++;
    } else {
      if(isnan(val))              node = casadi_limits<SX>::nan.node;
      else if(isinf(val))         node = val > 0 ? casadi_limits<SX>::inf.node : casadi_limits<SX>::minus_inf.node;
      else                        node = RealtypeSX::create(val);
      node->count++;
    }
  }

  SX::SX(const std::string& name){
    node = new SymbolicSX(name);  
    node->count++;
  }

  SX::~SX(){
    if(--node->count == 0) delete node;
  }

  SX& SX::operator=(const SX &scalar){
    // quick return if the old and new pointers point to the same object
    if(node == scalar.node) return *this;

    // decrease the counter and delete if this was the last pointer        
    if(--node->count == 0) delete node;

    // save the new pointer
    node = scalar.node;
    node->count++;
    return *this;
  }

  void SX::assignIfDuplicate(const SX& scalar, int depth){
    casadi_assert(depth>=1);
    if(!isEqual(scalar,0) && isEqual(scalar,depth)){
      *this = scalar;
    }
  }

  SXNode* SX::assignNoDelete(const SX& scalar){
    // Return value
    SXNode* ret = node;

    // quick return if the old and new pointers point to the same object
    if(node == scalar.node) return ret;

    // decrease the counter but do not delete if this was the last pointer
    --node->count;

    // save the new pointer
    node = scalar.node;
    node->count++;
  
    // Return a pointer to the old node
    return ret;
  }

  SX& SX::operator=(double scalar){
    return *this = SX(scalar);
  }

  std::ostream &operator<<(std::ostream &stream, const SX &scalar)
  {
    scalar.node->print(stream);  
    return stream;
  }

  void SX::print(std::ostream &stream, long& remaining_calls) const{
    if(remaining_calls>0){
      remaining_calls--;
      node->print(stream,remaining_calls);
    } else {
      stream << "...";
    }
  }

  SX SX::operator-() const{
    if(node->hasDep() && node->getOp() == OP_NEG)
      return node->dep(0);
    else if(node->isZero())
      return 0;
    else if(node->isMinusOne())
      return 1;
    else if(node->isOne())
      return -1;
    else
      return UnarySX::create(OP_NEG, *this);
  }

  SX SX::sign() const{
    return UnarySX::create(OP_SIGN, *this);
  }

  SX SX::erfinv() const{
    return UnarySX::create(OP_ERFINV,*this);
  }

  bool SX::__nonzero__() const {
    if (isConstant()) return !isZero();
    casadi_error("Cannot compute the truth value of a CasADi SX symbolic expression.")
      }

  SX SX::__add__(const SX& y) const{
    // NOTE: Only simplifications that do not result in extra nodes area allowed
    
    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_ADD,*this, y);
    
    if(node->isZero())
      return y;
    else if(y->isZero()) // term2 is zero
      return *this;
    else if(y.hasDep() && y.getOp()==OP_NEG) // x + (-y) -> x - y
      return __sub__(-y);
    else if(hasDep() && getOp()==OP_NEG) // (-x) + y -> y - x
      return y.__sub__(getDep());
    else if(hasDep() && getOp()==OP_MUL && 
            y.hasDep() && y.getOp()==OP_MUL && 
            getDep(0).isConstant() && getDep(0).getValue()==0.5 && 
            y.getDep(0).isConstant() && y.getDep(0).getValue()==0.5 &&
            y.getDep(1).isEqual(getDep(1),eq_depth_)) // 0.5x+0.5x = x
      return getDep(1);
    else if(hasDep() && getOp()==OP_DIV && 
            y.hasDep() && y.getOp()==OP_DIV && 
            getDep(1).isConstant() && getDep(1).getValue()==2 && 
            y.getDep(1).isConstant() && y.getDep(1).getValue()==2 &&
            y.getDep(0).isEqual(getDep(0),eq_depth_)) // x/2+x/2 = x
      return getDep(0);
    else if(hasDep() && getOp()==OP_SUB && getDep(1).isEqual(y,eq_depth_))
      return getDep(0);
    else if(y.hasDep() && y.getOp()==OP_SUB && isEqual(y.getDep(1),eq_depth_))
      return y.getDep(0);
    else // create a new branch
      return BinarySX::create(OP_ADD,*this, y);
  }

  SX SX::__sub__(const SX& y) const{
    // Only simplifications that do not result in extra nodes area allowed
    
    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_SUB,*this,y);
    
    if(y->isZero()) // term2 is zero
      return *this;
    if(node->isZero()) // term1 is zero
      return -y;
    if(isEqual(y,eq_depth_)) // the terms are equal
      return 0;
    else if(y.hasDep() && y.getOp()==OP_NEG) // x - (-y) -> x + y
      return __add__(-y);
    else if(hasDep() && getOp()==OP_ADD && getDep(1).isEqual(y,eq_depth_))
      return getDep(0);
    else if(hasDep() && getOp()==OP_ADD && getDep(0).isEqual(y,eq_depth_))
      return getDep(1);
    else if(y.hasDep() && y.getOp()==OP_ADD && isEqual(y.getDep(1),eq_depth_))
      return -y.getDep(0);
    else if(y.hasDep() && y.getOp()==OP_ADD && isEqual(y.getDep(0),eq_depth_))
      return -y.getDep(1);
    else // create a new branch
      return BinarySX::create(OP_SUB,*this,y);
  }

  SX SX::__mul__(const SX& y) const{
  
    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_MUL,*this,y);
    
    // Only simplifications that do not result in extra nodes area allowed
    if(y.isEqual(*this,eq_depth_))
      return sq();
    else if(!isConstant() && y.isConstant())
      return y.__mul__(*this);
    else if(node->isZero() || y->isZero()) // one of the terms is zero
      return 0;
    else if(node->isOne()) // term1 is one
      return y;
    else if(y->isOne()) // term2 is one
      return *this;
    else if(y->isMinusOne())
      return -(*this);
    else if(node->isMinusOne())
      return -y;
    else if(y.hasDep() && y.getOp()==OP_INV)
      return (*this)/y.inv();
    else if(hasDep() && getOp()==OP_INV)
      return y/inv();
    else if(isConstant() && y.hasDep() && y.getOp()==OP_MUL && y.getDep(0).isConstant() && getValue()*y.getDep(0).getValue()==1) // 5*(0.2*x) = x
      return y.getDep(1);
    else if(isConstant() && y.hasDep() && y.getOp()==OP_DIV && y.getDep(1).isConstant() && getValue()==y.getDep(1).getValue()) // 5*(x/5) = x
      return y.getDep(0);
    else if(hasDep() && getOp()==OP_DIV && getDep(1).isEqual(y,eq_depth_)) // ((2/x)*x)
      return getDep(0);
    else if(y.hasDep() && y.getOp()==OP_DIV && y.getDep(1).isEqual(*this,eq_depth_)) // ((2/x)*x)
      return y.getDep(0);
    else     // create a new branch
      return BinarySX::create(OP_MUL,*this,y);
  }


  bool SX::isDoubled() const{
    return isOp(OP_ADD) && node->dep(0).isEqual(node->dep(1),eq_depth_);
  }
    
  SX SX::__div__(const SX& y) const{
    // Only simplifications that do not result in extra nodes area allowed
    
    if (!CasadiOptions::simplification_on_the_fly) return BinarySX::create(OP_DIV,*this,y);

    if(y->isZero()) // term2 is zero
      return casadi_limits<SX>::nan;
    else if(node->isZero()) // term1 is zero
      return 0;
    else if(y->isOne()) // term2 is one
      return *this;
    else if(y->isMinusOne())
      return -(*this);
    else if(isEqual(y,eq_depth_)) // terms are equal
      return 1;
    else if(isDoubled() && y.isEqual(2))
      return node->dep(0);
    else if(isOp(OP_MUL) && y.isEqual(node->dep(0),eq_depth_))
      return node->dep(1);
    else if(isOp(OP_MUL) && y.isEqual(node->dep(1),eq_depth_))
      return node->dep(0);
    else if(node->isOne())
      return y.inv();
    else if(y.hasDep() && y.getOp()==OP_INV)
      return (*this)*y.inv();
    else if(isDoubled() && y.isDoubled())
      return node->dep(0) / y->dep(0);
    else if(y.isConstant() && hasDep() && getOp()==OP_DIV && getDep(1).isConstant() && y.getValue()*getDep(1).getValue()==1) // (x/5)/0.2 
      return getDep(0);
    else if(y.hasDep() && y.getOp()==OP_MUL && y.getDep(1).isEqual(*this,eq_depth_)) // x/(2*x) = 1/2
      return BinarySX::create(OP_DIV,1,y.getDep(0));
    else if(hasDep() && getOp()==OP_NEG && getDep(0).isEqual(y,eq_depth_))      // (-x)/x = -1
      return -1;
    else if(y.hasDep() && y.getOp()==OP_NEG && y.getDep(0).isEqual(*this,eq_depth_))      // x/(-x) = 1
      return -1;
    else if(y.hasDep() && y.getOp()==OP_NEG && hasDep() && getOp()==OP_NEG && getDep(0).isEqual(y.getDep(0),eq_depth_))      // (-x)/(-x) = 1
      return 1;
    else if(isOp(OP_DIV) && y.isEqual(node->dep(0),eq_depth_))
      return node->dep(1).inv();
    else // create a new branch
      return BinarySX::create(OP_DIV,*this,y);
  }

  SX SX::inv() const{
    if(node->hasDep() && node->getOp()==OP_INV){
      return node->dep(0);
    } else {
      return UnarySX::create(OP_INV,*this);
    }
  }

  Matrix<SX> SX::fmin(const Matrix<SX>& b) const { 
    return Matrix<SX>(*this).fmin(b);
  }
  Matrix<SX> SX::fmax(const Matrix<SX>& b) const {
    return Matrix<SX>(*this).fmax(b);
  }
  Matrix<SX> SX::constpow(const Matrix<SX>& n) const {
    return Matrix<SX>(*this).__constpow__(n);
  }

  Matrix<SX> SX::arctan2(const Matrix<SX>& b) const { 
    return Matrix<SX>(*this).arctan2(b);
  }

  SX SX::__le__(const SX& y) const{
    if((y-(*this)).isNonNegative())
      return 1;
    else
      return BinarySX::create(OP_LE,*this,y);
  }

  SX SX::__lt__(const SX& y) const{
    if(((*this)-y).isNonNegative())
      return 0;
    else
      return BinarySX::create(OP_LT,*this,y);
  }

  SX SX::__eq__(const SX& y) const{
    if(isEqual(y))
      return 1;
    else
      return BinarySX::create(OP_EQ,*this,y);
  }

  SX SX::__ne__(const SX& y) const{
    if(isEqual(y))
      return 0;
    else
      return BinarySX::create(OP_NE,*this,y);
  }

  SXNode* SX::get() const{
    return node;
  }

  const SXNode* SX::operator->() const{
    return node;
  }

  SXNode* SX::operator->(){
    return node;
  }

  SX if_else(const SX& cond, const SX& if_true, const SX& if_false){
    return if_else_zero(cond,if_true) + if_else_zero(!cond,if_false);
  }

  SX SX::binary(int op, const SX& x, const SX& y){
    return BinarySX::create(Operation(op),x,y);    
  }

  SX SX::unary(int op, const SX& x){
    return UnarySX::create(Operation(op),x);  
  }

  // SX::operator vector<SX>() const{
  //   vector<SX> ret(1);
  //   ret[0] = *this;
  //   return ret;
  // }

  string SX::toString() const{
    stringstream ss;
    ss << *this;
    return ss.str();
  }

  bool SX::isLeaf() const {
    if (!node) return true;
    return node->isConstant() || node->isSymbolic();
  }

  bool SX::isCommutative() const{
    if (!hasDep()) throw CasadiException("SX::isCommutative: must be binary");
    return operation_checker<CommChecker>(getOp());
  }

  bool SX::isConstant() const{
    return node->isConstant();
  }

  bool SX::isInteger() const{
    return node->isInteger();
  }

  bool SX::isSymbolic() const{
    return node->isSymbolic();
  }

  bool SX::hasDep() const{
    return node->hasDep();
  }

  bool SX::isZero() const{
    return node->isZero();
  }
  
  bool SX::isAlmostZero(double tol) const{
    return node->isAlmostZero(tol);
  }

  bool SX::isOne() const{
    return node->isOne();
  }

  bool SX::isMinusOne() const{
    return node->isMinusOne();
  }

  bool SX::isNan() const{
    return node->isNan();
  }

  bool SX::isInf() const{
    return node->isInf();
  }

  bool SX::isMinusInf() const{
    return node->isMinusInf();
  }

  const std::string& SX::getName() const{
    return node->getName();
  }

  int SX::getOp() const{
    return node->getOp();
  }

  bool SX::isOp(int op) const{
    return hasDep() && op==getOp();
  }

  bool SX::isEqual(const SX& ex, int depth) const{
    if(node==ex.get())
      return true;
    else if(depth>0)
      return node->isEqual(ex.get(),depth);
    else
      return false;
  }

  bool SX::isNonNegative() const{
    if(isConstant())
      return getValue()>=0;
    else if(isOp(OP_SQ) || isOp(OP_FABS))
      return true;
    else
      return false;
  }

  double SX::getValue() const{
    return node->getValue();
  }

  int SX::getIntValue() const{
    return node->getIntValue();
  }

  SX SX::getDep(int ch) const{
    casadi_assert(ch==0 || ch==1;)
      return node->dep(ch);
  }

  int SX::getNdeps() const {
    if (!hasDep()) throw CasadiException("SX::getNdeps: must be binary");
    return casadi_math<double>::ndeps(getOp());
  }

  long SX::__hash__() const {
    if (!node) return 0;
    return (long) node;
  }

  template<>
  bool __nonzero__<SX>(const SX& val) { return val.__nonzero__();} 

  const SX casadi_limits<SX>::zero(new ZeroSX(),false); // node corresponding to a constant 0
  const SX casadi_limits<SX>::one(new OneSX(),false); // node corresponding to a constant 1
  const SX casadi_limits<SX>::two(IntegerSX::create(2),false); // node corresponding to a constant 2
  const SX casadi_limits<SX>::minus_one(new MinusOneSX(),false); // node corresponding to a constant -1
  const SX casadi_limits<SX>::nan(new NanSX(),false);
  const SX casadi_limits<SX>::inf(new InfSX(),false);
  const SX casadi_limits<SX>::minus_inf(new MinusInfSX(),false);

  bool casadi_limits<SX>::isZero(const SX& val){ 
    return val.isZero();
  }
  
  bool casadi_limits<SX>::isAlmostZero(const SX& val, double tol){ 
    return val.isAlmostZero(tol);
  }

  bool casadi_limits<SX>::isOne(const SX& val){ 
    return val.isOne();
  }

  bool casadi_limits<SX>::isMinusOne(const SX& val){ 
    return val.isMinusOne();
  }

  bool casadi_limits<SX>::isConstant(const SX& val){
    return val.isConstant();
  }

  bool casadi_limits<SX>::isInteger(const SX& val){
    return val.isInteger();
  }

  bool casadi_limits<SX>::isInf(const SX& val){
    return val.isInf();
  }

  bool casadi_limits<SX>::isMinusInf(const SX& val){
    return val.isMinusInf();
  }

  bool casadi_limits<SX>::isNaN(const SX& val){
    return val.isNan();
  }

  SX SX::exp() const{
    return UnarySX::create(OP_EXP,*this);
  }

  SX SX::log() const{
    return UnarySX::create(OP_LOG,*this);
  }

  SX SX::log10() const{
    return log()*(1/std::log(10.));
  }

  SX SX::sqrt() const{
    if(isOp(OP_SQ))
      return node->dep(0).fabs();
    else
      return UnarySX::create(OP_SQRT,*this);
  }

  SX SX::sq() const{
    if(isOp(OP_SQRT))
      return node->dep(0);
    else
      return UnarySX::create(OP_SQ,*this);
  }

  SX SX::sin() const{
    return UnarySX::create(OP_SIN,*this);
  }

  SX SX::cos() const{
    return UnarySX::create(OP_COS,*this);
  }

  SX SX::tan() const{
    return UnarySX::create(OP_TAN,*this);
  }

  SX SX::arcsin() const{
    return UnarySX::create(OP_ASIN,*this);
  }

  SX SX::arccos() const{
    return UnarySX::create(OP_ACOS,*this);
  }

  SX SX::arctan() const{
    return UnarySX::create(OP_ATAN,*this);
  }

  SX SX::sinh() const{
    if(node->isZero())
      return 0;
    else
      return UnarySX::create(OP_SINH,*this);
  }

  SX SX::cosh() const{
    if(node->isZero())
      return 1;
    else
      return UnarySX::create(OP_COSH,*this);
  }

  SX SX::tanh() const{
    if(node->isZero())
      return 0;
    else
      return UnarySX::create(OP_TANH,*this);
  }

  SX SX::arctanh() const{
    if(node->isZero())
      return 0;
    else
      return UnarySX::create(OP_ATANH,*this);
  }

  SX SX::arccosh() const{
    if(node->isOne())
      return 0;
    else
      return UnarySX::create(OP_ACOSH,*this);
  }

  SX SX::arcsinh() const{
    if(node->isZero())
      return 0;
    else
      return UnarySX::create(OP_ASINH,*this);
  }

  SX SX::floor() const{
    return UnarySX::create(OP_FLOOR,*this);
  }

  SX SX::ceil() const{
    return UnarySX::create(OP_CEIL,*this);
  }

  SX SX::erf() const{
    return UnarySX::create(OP_ERF,*this);
  }

  SX SX::fabs() const{
    if(isOp(OP_FABS) || isOp(OP_SQ))
      return *this;
    else
      return UnarySX::create(OP_FABS,*this);
  }

  SX::operator Matrix<SX>() const{
    return Matrix<SX>(1,1,*this);
  }

  SX SX::fmin(const SX &b) const{
    return BinarySX::create(OP_FMIN,*this,b);
  }

  SX SX::fmax(const SX &b) const{
    return BinarySX::create(OP_FMAX,*this,b);
  }

  SX SX::arctan2(const SX &b) const{
    return BinarySX::create(OP_ATAN2,*this,b);
  }

  SX SX::printme(const SX &b) const{
    return BinarySX::create(OP_PRINTME,*this,b);
  }

  SX SX::__pow__(const SX& n) const{
    if(n->isConstant()) {
      if (n->isInteger()){
        int nn = n->getIntValue();
        if(nn == 0)
          return 1;
        else if(nn>100 || nn<-100) // maximum depth
          return BinarySX::create(OP_CONSTPOW,*this,nn);
        else if(nn<0) // negative power
          return 1/pow(*this,-nn);
        else if(nn%2 == 1) // odd power
          return *this*pow(*this,nn-1);
        else{ // even power
          SX rt = pow(*this,nn/2);
          return rt*rt;
        }
      } else if(n->getValue()==0.5){
        return sqrt();
      } else {
        return BinarySX::create(OP_CONSTPOW,*this,n);
      }
    } else {
      return BinarySX::create(OP_POW,*this,n);
    }
  }

  SX SX::__constpow__(const SX& n) const{
    return BinarySX::create(OP_CONSTPOW,*this,n);
  }

  SX SX::constpow(const SX& n) const{
    return BinarySX::create(OP_CONSTPOW,*this,n);
  }

  SX SX::logic_not() const{
    if(hasDep() && getOp() == OP_NOT){
      return getDep();
    } else {
      return UnarySX::create(OP_NOT, *this);
    }
  }

  SX SX::logic_and(const SX& y) const{
    return BinarySX::create(OP_AND,*this,y);
  }

  SX SX::logic_or(const SX& y) const{
    return BinarySX::create(OP_OR,*this,y);
  }

  SX SX::if_else_zero(const SX& y) const{
    if(y->isZero()){
      return y;
    } else if(isConstant()){
      if(getValue()!=0) return y;
      else              return 0;
    } else {
      return BinarySX::create(OP_IF_ELSE_ZERO,*this,y);
    }
  }

  int SX::getTemp() const{
    return (*this)->temp;
  }
    
  void SX::setTemp(int t){
    (*this)->temp = t;
  }

  bool SX::marked() const{
    return (*this)->marked();
  }
    
  void SX::mark(){
    (*this)->mark();
  }

  long SX::max_num_calls_in_print_ = 10000;

  void SX::setMaxNumCallsInPrint(long num){
    max_num_calls_in_print_ = num;
  }

  long SX::getMaxNumCallsInPrint(){
    return max_num_calls_in_print_;
  }

  int SX::eq_depth_ = 1;

  void SX::setEqualityCheckingDepth(int eq_depth){
    eq_depth_ = eq_depth;
  }

  int SX::getEqualityCheckingDepth(){
    return eq_depth_;
  }

} // namespace CasADi

using namespace CasADi;
namespace std{

  SX numeric_limits<SX>::infinity() throw(){
    return CasADi::casadi_limits<SX>::inf;
  }

  SX numeric_limits<SX>::quiet_NaN() throw(){
    return CasADi::casadi_limits<SX>::nan;
  }

  SX numeric_limits<SX>::min() throw(){
    return SX(numeric_limits<double>::min());
  }

  SX numeric_limits<SX>::max() throw(){
    return SX(numeric_limits<double>::max());
  }

  SX numeric_limits<SX>::epsilon() throw(){
    return SX(numeric_limits<double>::epsilon());
  }

  SX numeric_limits<SX>::round_error() throw(){
    return SX(numeric_limits<double>::round_error());
  }


} // namespace std

