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
#include <stack>
#include <cassert>
#include "../pre_c99_support.hpp"

using namespace std;
namespace CasADi{


SX::SX(){
  node = casadi_limits<SX>::nan.node;
  node->count++;
}

SX::SX(SXNode* node_) : node(node_){
  node->count++;
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
    else                        node = new IntegerSXNode(intval);
    node->count++;
  } else {
    if(isnan(val))              node = casadi_limits<SX>::nan.node;
    else if(isinf(val))         node = val > 0 ? casadi_limits<SX>::inf.node : casadi_limits<SX>::minus_inf.node;
    else                        node = new RealtypeSXNode(val);
    node->count++;
  }
}

SX::SX(const std::string& name){
  node = new SymbolicSXNode(name);  
  node->count++;
}

SX::SX(const char name[]){
  node = new SymbolicSXNode(name);  
  node->count++;
}

SX::SX(const Matrix<SX>& m){
  // Assign the node a temporary value
  node = casadi_limits<SX>::zero.node;
  node->count++;
  
  // Copy the value
  *this = m.toScalar();
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

SX& operator+=(SX &ex, const SX &el){
  return ex = ex + el;
}

SX& operator-=(SX &ex, const SX &el){
  return ex = ex - el;
}

SX SX::operator-() const{
  if(node->hasDep() && node->getOp() == NEG)
    return node->dep(0);
  else if(node->isZero())
    return 0;
  else if(node->isMinusOne())
    return 1;
  else if(node->isOne())
    return -1;
  else
   return SX(new BinarySXNode(NEG, *this));
}

SX& operator*=(SX &ex, const SX &el){
 return ex = ex * el;
}

SX& operator/=(SX &ex, const SX &el){
  return ex = ex / el;
}

SX sign(const SX& x){
  return timesTwo(x>=0)-1;
}


SX SX::add(const SX& y) const{
  if(node->isZero())
    return y;
  else if(y->isZero()) // term2 is zero
    return *this;
  else if(y.isBinary() && y.getOp()==NEG) // x + (-y) -> x - y
    return sub(-y);
  else if(isBinary() && getOp()==NEG) // (-x) + y -> y - x
    return y.sub(getDep());
  else // create a new branch
    return SX(new BinarySXNode(ADD, *this, y));
}

SX SX::sub(const SX& y) const{
  if(y->isZero()) // term2 is zero
    return *this;
  if(node->isZero()) // term1 is zero
    return -y;
  if(node->isEqual(y)) // the terms are equal
    return 0;
  else if(y.isBinary() && y.getOp()==NEG) // x - (-y) -> x + y
    return add(-y);
  else // create a new branch
    return SX(new BinarySXNode(SUB, *this, y));
}

SX SX::mul(const SX& y) const{
  if(node->isZero() || y->isZero()) // one of the terms is zero
    return 0;
  else if(node->isOne()) // term1 is one
    return y;
  else if(y->isOne()) // term2 is one
    return *this;
  else if(y->isMinusOne())
    return -(*this);
  else if(node->isMinusOne())
    return -y;
  else if(y.isBinary() && y.getOp()==INV)
    return (*this)/y.inv();
  else if(isBinary() && getOp()==INV)
    return y/inv();
  else     // create a new branch
    return SX(new BinarySXNode(MUL,*this,y));
}

SX SX::div(const SX& y) const{
  if(y->isZero()) // term2 is zero
    return casadi_limits<SX>::nan;
  else if(node->isZero()) // term1 is zero
    return 0;
  else if(y->isOne()) // term2 is one
    return *this;
  else if(node->isEqual(y)) // terms are equal
    return 1;
  else if(y.isEqual(2) && node->hasDep() && node->getOp()==ADD && node->dep(0).isEqual(node->dep(1)))
    return node->dep(0);
  else if(node->isOne())
    return y.inv();
  else if(y.isBinary() && y.getOp()==INV)
    return (*this)*y.inv();
  else // create a new branch
    return SX(new BinarySXNode(DIV,*this,y));
}

SX SX::inv() const{
  if(node->hasDep() && node->getOp()==INV){
    return node->dep(0);
  } else {
    return SX(new BinarySXNode(INV,*this));
  }
}

SX operator+(const SX &x, const SX &y){
  return x.add(y);
}

SX operator-(const SX &x, const SX &y){
  return x.sub(y);
}

SX operator*(const SX &x, const SX &y){
  return x.mul(y);
}


SX operator/(const SX &x, const SX &y){
  return x.div(y);
}

SX operator<=(const SX &a, const SX &b){
  return b>=a;
}

SX operator>=(const SX &a, const SX &b){
  // Move everything to one side
  SX x = a-b;
  if(x->isConstant())
    return x->getValue()>=0; // ok since the result will be either 0 or 1, i.e. no new nodes
  else
    return SX(new BinarySXNode(STEP,x));
}

SX operator<(const SX &a, const SX &b){
  return !(a>=b);
}

SX operator>(const SX &a, const SX &b){
  return !(a<=b);
}

SX operator&&(const SX &a, const SX &b){
  return a+b>=2;
}

SX operator||(const SX &a, const SX &b){
  return !(!a && !b);
}

SX operator==(const SX &x, const SX &y){
  if(x.isEqual(y))
    return 1; // this also covers the case when both are constant and equal
  else if(x.isConstant() && y.isConstant())
    return 0;
  else // create a new node
    return SX(new BinarySXNode(EQUALITY,x,y));
}

SX operator!=(const SX &a, const SX &b){
  return !(a == b);
}

SX operator!(const SX &a){
  return 1-a;
}

SXNode* const SX::get() const{
  return node;
}

const SXNode* SX::operator->() const{
  return node;
}

SXNode* SX::operator->(){
  return node;
}

SX if_else(const SX& cond, const SX& if_true, const SX& if_false){
  return if_false + (if_true-if_false)*cond;
}

SX SX::binary(int op, const SX& x, const SX& y){
  return SX(new BinarySXNode(Operation(op),x,y));    
}

SX SX::unary(int op, const SX& x){
  return SX(new BinarySXNode(Operation(op),x));  
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

bool SX::isConstant() const{
  return node->isConstant();
}

bool SX::isInteger() const{
  return node->isInteger();
}

bool SX::isSymbolic() const{
  return node->isSymbolic();
}

bool SX::isBinary() const{
  return node->hasDep();
}

bool SX::isZero() const{
  return node->isZero();
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

bool SX::isEqual(const SX& scalar) const{
  return node->isEqual(scalar);
}

double SX::getValue() const{
  return node->getValue();
}

int SX::getIntValue() const{
  return node->getIntValue();
}

SX SX::getDep(int ch) const{
  return node->dep(ch);
}

const SX casadi_limits<SX>::zero(new ZeroSXNode()); // node corresponding to a constant 0
const SX casadi_limits<SX>::one(new OneSXNode()); // node corresponding to a constant 1
const SX casadi_limits<SX>::two(new IntegerSXNode(2)); // node corresponding to a constant 2
const SX casadi_limits<SX>::minus_one(new MinusOneSXNode()); // node corresponding to a constant -1
const SX casadi_limits<SX>::nan(new NanSXNode());
const SX casadi_limits<SX>::inf(new InfSXNode());
const SX casadi_limits<SX>::minus_inf(new MinusInfSXNode());

bool casadi_limits<SX>::isZero(const SX& val){ 
  return val.isZero();
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
  return SX(new BinarySXNode(EXP,*this));
}

SX SX::log() const{
  return SX(new BinarySXNode(LOG,*this));
}

SX SX::log10() const{
  return log()*(1/std::log(10.));
}

SX SX::sqrt() const{
  if(isOne() || isZero())
    return *this;
  else if(isBinary() && getOp()==MUL && getDep(0).isEqual(getDep(1)))
    return getDep().fabs();
  else
    return SX(new BinarySXNode(SQRT,*this));
}

SX SX::sin() const{
  if(node->isZero())
    return 0;
  else
    return SX(new BinarySXNode(SIN,*this));
}

SX SX::cos() const{
  if(node->isZero())
    return 1;
  else
    return SX(new BinarySXNode(COS,*this));
}

SX SX::tan() const{
  if(node->isZero())
    return 0;
  else
    return SX(new BinarySXNode(TAN,*this));
}

SX SX::arcsin() const{
  return SX(new BinarySXNode(ASIN,*this));
}

SX SX::arccos() const{
  return SX(new BinarySXNode(ACOS,*this));
}

SX SX::arctan() const{
  return SX(new BinarySXNode(ATAN,*this));
}

SX SX::sinh() const{
  if(node->isZero())
    return 0;
  else
    return SX(new BinarySXNode(SINH,*this));
}

SX SX::cosh() const{
  if(node->isZero())
    return 1;
  else
    return SX(new BinarySXNode(COSH,*this));
}

SX SX::tanh() const{
  if(node->isZero())
    return 0;
  else
    return SX(new BinarySXNode(TANH,*this));
}

SX SX::floor() const{
  return SX(new BinarySXNode(FLOOR,*this));
}

SX SX::ceil() const{
  return SX(new BinarySXNode(CEIL,*this));
}

SX SX::erf() const{
  return SX(new BinarySXNode(ERF,*this));
}

SX SX::fabs() const{
  if(isConstant() && getValue()>=0)
    return *this;
  else
    return sign(*this)**this;
}

SX casadi_operators<SX>::add(const SX&x, const SX&y){
  return x.add(y);
}

SX casadi_operators<SX>::sub(const SX&x, const SX&y){
  return x.sub(y);
}

SX casadi_operators<SX>::mul(const SX&x, const SX&y){
  return x.mul(y);
}

SX casadi_operators<SX>::div(const SX&x, const SX&y){
  return x.div(y);
}

SX casadi_operators<SX>::fmin(const SX&x, const SX&y){
  return x.fmin(y);
}

SX casadi_operators<SX>::fmax(const SX&x, const SX&y){
  return x.fmax(y);
}

SX casadi_operators<SX>::pow(const SX&x, const SX&y){
  return x.__pow__(y);
}

SX casadi_operators<SX>::constpow(const SX&x, const SX&y){
  return x.__constpow__(y);
}

SX casadi_operators<SX>::printme(const SX&x, const SX&y){
  return x.printme(y);
}

SX casadi_operators<SX>::sin(const SX&x){ 
  return x.sin();
}

SX casadi_operators<SX>::cos(const SX&x){ 
  return x.cos();
}

SX casadi_operators<SX>::tan(const SX&x){ 
  return x.tan();
}

SX casadi_operators<SX>::asin(const SX&x){ 
  return x.arcsin();
}

SX casadi_operators<SX>::acos(const SX&x){ 
  return x.arccos();
}

SX casadi_operators<SX>::atan(const SX&x){ 
  return x.arctan();
}

SX casadi_operators<SX>::sinh(const SX&x){ 
  return x.sinh();
}

SX casadi_operators<SX>::cosh(const SX&x){ 
  return x.cosh();
}

SX casadi_operators<SX>::tanh(const SX&x){ 
  return x.tanh();
}

SX casadi_operators<SX>::neg(const SX&x){
  return x.operator-();
}

SX casadi_operators<SX>::log(const SX&x){ 
  return x.log();
}

SX casadi_operators<SX>::exp(const SX&x){ 
  return x.exp();
}

SX casadi_operators<SX>::sqrt(const SX&x){ 
  return x.sqrt();
}

SX casadi_operators<SX>::floor(const SX&x){ 
  return x.floor();
}

SX casadi_operators<SX>::ceil(const SX&x){ 
  return x.ceil();
}

SX casadi_operators<SX>::fabs(const SX&x){ 
  return x.fabs();
}

SX casadi_operators<SX>::erf(const SX&x){ 
  return x.erf();
}

SX::operator Matrix<SX>() const{
  return Matrix<SX>(1,1,*this);
}

SX SX::fmin(const SX &b) const{
  return SX(new BinarySXNode(FMIN,*this,b));
}

SX SX::fmax(const SX &b) const{
  return SX(new BinarySXNode(FMAX,*this,b));
}

SX SX::printme(const SX &b) const{
  return SX(new BinarySXNode(PRINTME,*this,b));
}

SX SX::__pow__(const SX& n) const{
  if(n->isConstant()) {
    if (n->isInteger()){
      int nn = n->getIntValue();
      if(nn == 0)
        return 1;
      else if(nn>100 || nn<-100) // maximum depth
        return SX(new BinarySXNode(CONSTPOW,*this,nn));
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
      return SX(new BinarySXNode(CONSTPOW,*this,n));
    }
  } else {
    return SX(new BinarySXNode(POW,*this,n));
  }
}

SX SX::__constpow__(const SX& n) const{
  return SX(new BinarySXNode(CONSTPOW,*this,n));
}

SX SX::constpow(const SX& n) const{
  return SX(new BinarySXNode(CONSTPOW,*this,n));
}

int SX::getTemp() const{
  return (*this)->temp;
}
    
void SX::setTemp(int t){
  (*this)->temp = t;
}

long SX::max_num_calls_in_print = 10000;


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

