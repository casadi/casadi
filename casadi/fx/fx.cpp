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

#include "fx.hpp"
#include "../mx/evaluation.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"

using namespace std;

namespace CasADi{
  
FX::FX(){
}

FX::~FX(){
}

FX::FX(const FX& fx) : OptionsFunctionality(fx){
}
    
const FXNode* FX::operator->() const{
  return (const FXNode*)OptionsFunctionality::operator->();
}

FXNode* FX::operator->(){
  return (FXNode*)OptionsFunctionality::operator->();
}

MX FX::operator()(const MX &x, int ind) const{
  vector<MX> dep(1);
  dep[0] = x;
  MX ret;
  ret.assignNode(new Evaluation(*this,dep,ind));
  return ret;
}

MX FX::operator()(const vector<MX> &x, int ind) const{
  MX ret;
  ret.assignNode(new Evaluation(*this,x,ind));
  return ret;
}

void FX::evaluate(int fsens_order, int asens_order){
  (*this)->evaluate(fsens_order,asens_order);
}

void FX::solve(){
  (*this)->evaluate(0,0);
}

int FX::getNumInputs() const{
  return (*this)->input_.size();
}

int FX::getNumOutputs() const{
  return (*this)->output_.size();  
}

void FX::setNumInputs(int num_in){
  return (*this)->input_.resize(num_in);
}

void FX::setNumOutputs(int num_out){
  return (*this)->output_.resize(num_out);  
}


FX FX::jacobian(int iind, int oind){
  return (*this)->jacobian(iind,oind);  
}

FX FX::hessian(int iind, int oind){
  return (*this)->hessian(iind,oind);  
}

void FX::assertNode() const{
  if(!dynamic_cast<const FXNode*>(get()))
    throw CasadiException("FX::assertNode");
}




FXNode::FXNode(){
  setOption("name",            "unnamed_function"); // name of the function
  addOption("ad_order",          OT_INTEGER,   0); // use sparse jacobian
  addOption("sparse",            OT_BOOLEAN,   true); // function is sparse
  addOption("number_of_fwd_dir", OT_INTEGER,  1); // number of forward derivatives
  addOption("number_of_adj_dir", OT_INTEGER,  1); // number of adjoint derivatives

  ad_order_ = 0;
  is_init_ = false;
}

FXNode::~FXNode(){
}

void FXNode::init(){
  ad_order_ = getOption("ad_order").toInt();
  nfdir_ = getOption("number_of_fwd_dir").toInt();
  nadir_ = getOption("number_of_adj_dir").toInt();

  for(vector<FunctionIO>::iterator it=input_.begin(); it!=input_.end(); ++it){
    if(ad_order_==1){
      it->setNumFwdDir(nfdir_);
      it->setNumAdjDir(nadir_);
    }
    it->init();
  }

  for(vector<FunctionIO>::iterator it=output_.begin(); it!=output_.end(); ++it){
    if(ad_order_==1){
      it->setNumFwdDir(nfdir_);
      it->setNumAdjDir(nadir_);
    }
    it->init();
  }

  is_init_ = true;
}

void FXNode::print(ostream &stream) const{
  stream << "function(\"" << getOption("name") << "\")";
}

FX FXNode::jacobian(int iind, int oind){
  FX fcn;
  fcn.assignNode(this);
  return Jacobian(fcn,iind,oind);
}

FX FXNode::hessian(int iind, int oind){
  stringstream ss;
  ss << "FXNode::hessian: hessian not defined for class " << typeid(*this).name();
  throw CasadiException(ss.str());
}

void FXNode::assertInit() const{
  if(!is_init_)
    throw CasadiException("FXNode::assertInit: function has not been initialized");
}


FunctionIO& FX::input(int i){
  return (*this)->input(i);
}

const FunctionIO& FX::input(int i) const{
  return (*this)->input(i);
}
  
FunctionIO& FX::output(int i){
  return (*this)->output(i);
}

const FunctionIO& FX::output(int i) const{
  return (*this)->output(i);
}

FunctionIO& FXNode::input(int i){
  if(i>=input_.size()) throw CasadiException("FXNode::input: out of bounds");
  return input_.at(i);
}

const FunctionIO& FXNode::input(int i) const{
  if(i>=input_.size()) throw CasadiException("FXNode::input: out of bounds");
  return input_.at(i);  
}
  
FunctionIO& FXNode::output(int i){
  if(i>=output_.size()) throw CasadiException("FXNode::output: out of bounds");
  return output_.at(i);  
}

const FunctionIO& FXNode::output(int i) const{
  if(i>=output_.size()) throw CasadiException("FXNode::output: out of bounds");
  return output_.at(i);    
}




// void setv(double val, vector<double>& v){
//   if(v.size() != 1) throw CasadiException("setv(double,vector<double>&): dimension mismatch");
//   v[0] = val;
// }
// 
// void setv(int val, vector<double>& v){
//   setv(double(val),v);
// }
// 
// void setv(const double* val, vector<double>& v){
//   // no checking
//   copy(val,val+v.size(),v.begin());
// }
// 
// void setv(const vector<double>& val, vector<double>& v){
//   if(v.size() != val.size()) throw CasadiException("setv(const vector<double>&,vector<double>&): dimension mismatch");
//   copy(val.begin(),val.end(),v.begin());
// }
//   
// void getv(double &val, const vector<double>& v){
//   if(v.size() != 1) throw CasadiException("getv(double&,vector<double>&): dimension mismatch");
//   val = v[0];  
// }
// 
// void getv(double* val, const vector<double>& v){
//   // no checking
//   copy(v.begin(),v.end(),val);
// }
// 
// void getv(vector<double>& val, const vector<double>& v){
//   if(v.size() != val.size()) throw CasadiException("getv(vector<double>&,const vector<double>&): dimension mismatch");
//   copy(v.begin(),v.end(),val.begin());
// }

} // namespace CasADi

