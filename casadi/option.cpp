/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "option.hpp"
#include "stl_vector_tools.hpp"
#include "casadi_exception.hpp"

using namespace std;

namespace CasADi{


Option::Option(){
}

bool OptionNode::toBool() const{
  return bool(toInt());
}

int OptionNode::toInt() const{
  return toIntVector().at(0);
}

double OptionNode::toDouble() const{
  return toDoubleVector().at(0);
}

const vector<int>& OptionNode::toIntVector() const{
  if(n != i_vec.size()) throw CasadiException("OptionNode::toIntVector");
  return i_vec;
}

const vector<double>& OptionNode::toDoubleVector() const{
  if(n != d_vec.size()) throw CasadiException("OptionNode::toDoubleVector");
  return d_vec;
}

const string& OptionNode::toString() const{
  return str;
}

OptionNode::OptionNode(const vector<int>& i_vec_) : i_vec(i_vec_){
  stringstream ss;
  ss << i_vec;
  str = ss.str();
  is_string = true;
  n = i_vec.size();
  d_vec.resize(n);
  copy(i_vec.begin(),i_vec.end(),d_vec.begin());
}

OptionNode::OptionNode(const vector<double>& d_vec_) : d_vec(d_vec_){
  stringstream ss;
  ss << i_vec;
  str = ss.str();
  is_string = true;
  n = d_vec.size();
  i_vec.resize(n);
  copy(d_vec.begin(),d_vec.end(),i_vec.begin());
}

OptionNode::OptionNode(const string& s) : str(s){
  is_string = true;
  n = 1;
}


ostream& operator<<(ostream &stream, const Option& ref){
  if(!ref->d_vec.empty())
    return stream << ref->d_vec;
  else
    return stream << ref->str;
}

Option::Option(int i){
  vector<int> v(1);
  v[0] = i;
  assignNode(new OptionNode(v));
}

Option::Option(double d){
  vector<double> v(1);
  v[0] = d;
  assignNode(new OptionNode(v));
}

Option::Option(const vector<int>& iv){
  assignNode(new OptionNode(iv));
}

Option::Option(const vector<bool>& b_vec){
  vector<int> i_vec(b_vec.size());
  copy(b_vec.begin(),b_vec.end(), i_vec.begin());
  assignNode(new OptionNode(i_vec));
}

Option::Option(const vector<double>& dv){
  assignNode(new OptionNode(dv));
}

Option::Option(const string& s){
  assignNode(new OptionNode(s));
}

Option::Option(const char s[]){ 
  assignNode(new OptionNode(s));
}


OptionNode* Option::operator->(){
  return (OptionNode*)SharedObject::operator->();
}

const OptionNode* Option::operator->() const{
  return (const OptionNode*)SharedObject::operator->();
}

bool Option::toBool() const{
  return (*this)->toBool();
}

int Option::toInt() const{
  return (*this)->toInt();  
}

double Option::toDouble() const{
  return (*this)->toDouble();    
}

string Option::toString() const{
  return (*this)->toString();    
}

const vector<int>& Option::toIntVector() const{
  return (*this)->toIntVector();    
}

const vector<double>& Option::toDoubleVector() const{
  return (*this)->toDoubleVector();      
}


bool operator==(const Option& op1, const Option& op2){
  return !(op1 != op2);
}

bool operator!=(const Option& op1, const Option& op2){
  if(op1->is_string){
    if(!op2->is_string) return true;
    return op1.toString().compare(op2.toString()) != 0;
  } else {
    if(op2->is_string) return true;
    const vector<double> &v1 = op1.toDoubleVector();
    const vector<double> &v2 = op2.toDoubleVector();
    if(v1.size() != v2.size()) return true;
    for(int i=0; i<v1.size(); ++i)
      if(v1[i] != v2[i]) return true;
  }
  
  return false;
}



} // namespace CasADi

