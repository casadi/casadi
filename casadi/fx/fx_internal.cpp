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

#include "fx_internal.hpp"
#include "../mx/evaluation.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"

using namespace std;

namespace CasADi{
  
FXInternal::FXInternal(){
  setOption("name",            "unnamed_function"); // name of the function
  addOption("sparse",            OT_BOOLEAN,   true); // function is sparse
  addOption("number_of_fwd_dir", OT_INTEGER,  1); // number of forward derivatives
  addOption("number_of_adj_dir", OT_INTEGER,  1); // number of adjoint derivatives
  addOption("verbose",           OT_BOOLEAN,   false); // verbose evaluation -- for debugging
  is_init_ = false;
  verbose_ = false;
}

FXInternal::~FXInternal(){
}

void FXInternal::init(){
  nfdir_ = getOption("number_of_fwd_dir").toInt();
  nadir_ = getOption("number_of_adj_dir").toInt();
  verbose_ = getOption("verbose").toInt();

  for(vector<FunctionIO>::iterator it=input_.begin(); it!=input_.end(); ++it){
    it->dataF.resize(nfdir_);
    it->dataA.resize(nadir_);
    it->init();
  }

  for(vector<FunctionIO>::iterator it=output_.begin(); it!=output_.end(); ++it){
    it->dataF.resize(nfdir_);
    it->dataA.resize(nadir_);
    it->init();
  }

  is_init_ = true;
}

void FXInternal::print(ostream &stream) const{
  stream << "function(\"" << getOption("name") << "\")";
}

FX FXInternal::jacobian(int iind, int oind){
  FX fcn;
  fcn.assignNode(this);
  return Jacobian(fcn,iind,oind);
}

FX FXInternal::hessian(int iind, int oind){
  stringstream ss;
  ss << "FXInternal::hessian: hessian not defined for class " << typeid(*this).name();
  throw CasadiException(ss.str());
}

bool FXInternal::isInit() const{
  return is_init_;
}

FunctionIO& FXInternal::inputStruct(int i){
  if(i<0 || i>=input_.size()){
    stringstream ss;
    ss << "In function " << getOption("name") << ": input " << i << " not in interval [0," << input_.size() << "]"; 
    throw CasadiException(ss.str());
  }
  return input_.at(i);
}

const FunctionIO& FXInternal::inputStruct(int i) const{
  return const_cast<FXInternal*>(this)->inputStruct(i);
}
  
FunctionIO& FXInternal::outputStruct(int i){
  if(i<0 || i>=output_.size()){
    stringstream ss;
    ss << "In function " << getOption("name") << ": output " << i << " not in interval [0," << output_.size() << "]"; 
    throw CasadiException(ss.str());
  }
  return output_.at(i);
}

const FunctionIO& FXInternal::outputStruct(int i) const{
  return const_cast<FXInternal*>(this)->outputStruct(i);
}

void FXInternal::log(const std::string& msg) const{
  if(verbose()){
    cout << "CasADi log message: " << msg << endl;
  }
}

void FXInternal::log(const std::string& fcn, const std::string& msg) const{
  if(verbose()){
    cout << "CasADi log message: In \"" << fcn << "\" --- " << msg << endl;
  }
}

bool FXInternal::verbose() const{
  return verbose_;
}

bool FXInternal::monitored(const std::string& mod) const{
  return monitors_.count(mod)>0;
}

Matrix<double>& FXInternal::input(int iind){
  return inputStruct(iind).data;
}
    
const Matrix<double>& FXInternal::input(int iind) const{
  return inputStruct(iind).data;
}

Matrix<double>& FXInternal::output(int oind){
  return outputStruct(oind).data;
}
    
const Matrix<double>& FXInternal::output(int oind) const{
  return outputStruct(oind).data;
}

Matrix<double>& FXInternal::fwdSeed(int iind, int dir){
  return inputStruct(iind).dataF.at(dir);
}
    
const Matrix<double>& FXInternal::fwdSeed(int iind, int dir) const{
  return inputStruct(iind).dataF.at(dir);
}

Matrix<double>& FXInternal::fwdSens(int oind, int dir){
  return outputStruct(oind).dataF.at(dir);
}
    
const Matrix<double>& FXInternal::fwdSens(int oind, int dir) const{
  return outputStruct(oind).dataF.at(dir);
}

Matrix<double>& FXInternal::adjSeed(int oind, int dir){
  return outputStruct(oind).dataA.at(dir);
}
    
const Matrix<double>& FXInternal::adjSeed(int oind, int dir) const{
  return outputStruct(oind).dataA.at(dir);
}

Matrix<double>& FXInternal::adjSens(int iind, int dir){
  return inputStruct(iind).dataA.at(dir);
}
    
const Matrix<double>& FXInternal::adjSens(int iind, int dir) const{
  return inputStruct(iind).dataA.at(dir);
}

void FXInternal::setNumInputs(int num_in){
  return input_.resize(num_in);
}

void FXInternal::setNumOutputs(int num_out){
  return output_.resize(num_out);  
}

int FXInternal::getNumInputs() const{
  return input_.size();
}

int FXInternal::getNumOutputs() const{
  return output_.size();
}

const Dictionary & FXInternal::getStats() const {
  return stats_;
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

