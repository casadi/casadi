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

namespace CasADi{

const SX SX::zero(new ZeroSXNode()); // node corresponding to a constant 0
const SX SX::one(new OneSXNode()); // node corresponding to a constant 1
const SX SX::two(new IntegerSXNode(2)); // node corresponding to a constant 2
const SX SX::mone(new MinusOneSXNode()); // node corresponding to a constant -1
const SX SX::nan(new NanSXNode());
const SX SX::inf(new InfSXNode());
const SX SX::minf(new MinusInfSXNode());

SX::SX(){
  node = nan.node;
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
    if(intval == 0)             node = zero.node;
    else if(intval == 1)        node = one.node;
    else if(intval == 2)        node = two.node;
    else if(intval == -1)       node = mone.node;
    else                        node = new IntegerSXNode(intval);
    node->count++;
  } else {	
    if(ISNAN(val))              node = nan.node;
    else if(ISINF(val))         node = val > 0 ? inf.node : minf.node;
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

SX& operator+=(SX &ex, const SX &el){
  return ex = ex + el;
}

SX& operator-=(SX &ex, const SX &el){
  return ex = ex - el;
}

SX operator-(const SX &ex){
  if(ex->isBinary() && ex->getOp() == NEG_NODE)
    return ex->dependent(0);
  else if(ex->isMinusOne())
    return 1;
  else if(ex->isOne())
    return -1;
  else
   return SX(new BinarySXNode(NEG_NODE, ex));
}

SX& operator*=(SX &ex, const SX &el){
 return ex = ex * el;
}

SX& operator/=(SX &ex, const SX &el){
  return ex = ex / el;
}

SX sign(const SX& x){
  return 2*(x>=0)-1;
}


SX operator+(const SX &x, const SX &y){
    if(x->isZero())
      return y;
    else if(y->isZero()) // term2 is zero
      return x;
    else // create a new branch
      return SX(new BinarySXNode(ADD_NODE, x, y));
}

SX operator-(const SX &x, const SX &y){
    if(y->isZero()) // term2 is zero
      return x;
    if(x->isZero()) // term1 is zero
      return -y;
    if(x->isEqual(y)) // the terms are equal
      return 0;
    else // create a new branch
      return SX(new BinarySXNode(SUB_NODE, x, y));
}

SX operator*(const SX &x, const SX &y){
  if(x->isZero() || y->isZero()) // one of the terms is zero
      return 0;
    else if(x->isOne()) // term1 is one
      return y;
    else if(y->isOne()) // term2 is one
      return x;
    else if(y->isMinusOne())
      return -x;
    else if(x->isMinusOne())
      return -y;
    else     // create a new branch
      return SX(new BinarySXNode(MUL_NODE,x,y));
}


SX operator/(const SX &x, const SX &y){
    if(y->isZero()) // term2 is zero
      return SX::nan;
    else if(x->isZero()) // term1 is zero
      return 0;
    else if(y->isOne()) // term2 is one
      return x;
    else if(x->isEqual(y)) // terms are equal
      return 1;
    else // create a new branch
      return SX(new BinarySXNode(DIV_NODE,x,y));
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
    return SX(new BinarySXNode(STEP_NODE,x));
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
    if(x->isConstant())
    return x==y; // ok since the result will be either 0 or 1, i.e. no new nodes
  else
    return SX(new BinarySXNode(EQUALITY_NODE,x,y));
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
  return SX(new BinarySXNode(OPERATION(op),x,y));    
}

SX SX::unary(int op, const SX& x){
  return SX(new BinarySXNode(OPERATION(op),x));  
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

void make_symbolic(SX& v, const std::string& name){
  v = SX(name);
}

vector<SX> create_symbolic(const string& name, int n){
  vector<SX> ret(n);
  make_symbolic(ret,name);
  return ret;
}

vector< vector<SX> > create_symbolic(const std::string& name, int n, int m){
  vector< vector<SX> > ret(n,vector<SX>(m));
  make_symbolic(ret,name);
  return ret;
}

vector< vector< vector<SX> > > create_symbolic(const std::string& name, int n, int m, int p){
  vector< vector< vector<SX> > > ret(n,vector< vector<SX> >(m, vector<SX>(p)));
  make_symbolic(ret,name);
  return ret;
}

} // namespace CasADi

using namespace CasADi;
namespace std{

SX fabs(const SX& x){
  return sign(x)*x;
}

SX exp(const SX& x){
  return SX(new BinarySXNode(EXP_NODE,x));
}

SX log(const SX& x){
  return SX(new BinarySXNode(LOG_NODE,x));
}

SX sqrt(const SX& x){
  return SX(new BinarySXNode(SQRT_NODE,x));
}

SX sin(const SX& x){
  return SX(new BinarySXNode(SIN_NODE,x));
}

SX cos(const SX& x){;  
  return SX(new BinarySXNode(COS_NODE,x));
}

SX tan(const SX& x){
  return SX(new BinarySXNode(TAN_NODE,x));
}

SX atan(const SX& x){
  return SX(new BinarySXNode(ATAN_NODE,x));
}

SX asin(const SX& x){
  return SX(new BinarySXNode(ASIN_NODE,x));
}

SX acos(const SX& x){
  return SX(new BinarySXNode(ACOS_NODE,x));
}


SX pow(const SX& x, const SX& n){
  if (n->isInteger()){
    int nn = n->getIntValue();
    if(nn == 0)
      return 1;
    else if(nn<0) // negative power
      return 1/pow(x,-nn);
    else if(nn%2 == 1) // odd power
      return x*pow(x,nn-1);
    else{ // even power
      SX rt = pow(x,nn/2);
      return rt*rt;
    }
  } if(n->isConstant() && n->getValue()==0.5) {
    return sqrt(x);
  } else {
    return SX(new BinarySXNode(POW_NODE,x,n));
  }
}

SX erf(const SX& x){
  return SX(new BinarySXNode(ERF_NODE,x));
}

SX abs(const SX& x){
  return sign(x)*x;
}

SX fmin(const SX &a, const SX &b){
  return SX(new BinarySXNode(FMIN_NODE,a,b));
}

SX fmax(const SX &a, const SX &b){
  return SX(new BinarySXNode(FMAX_NODE,a,b));
}

SX floor(const SX& x){
  return SX(new BinarySXNode(FLOOR_NODE,x));
}

SX ceil(const SX& x){
  return SX(new BinarySXNode(CEIL_NODE,x));
}



} // namespace std

