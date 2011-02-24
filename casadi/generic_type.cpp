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

#include "generic_type_internal.hpp"
#include "stl_vector_tools.hpp"
#include "casadi_exception.hpp"
#include <cmath>

using namespace std;

namespace CasADi{

  typedef GenericTypeInternal<std::string> StringType;
  typedef GenericTypeInternal<double> DoubleType;
  typedef GenericTypeInternal<int> IntType;
  typedef GenericTypeInternal<std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<std::vector<int> > IntVectorType;

bool GenericType::isBool() const{
  return isInt() && (toInt()==0 || toInt()==1);
}

bool GenericType::isInt() const{
  return is_a<int>();
}
    
bool GenericType::isDouble() const{
  return is_a<double>();
}
    
bool GenericType::isString() const{
  return is_a<string>();
}

bool GenericType::isIntVector() const{
  return is_a<vector<int> >();
}
    
bool GenericType::isDoubleVector() const{
  return is_a<vector<double> >();
}
  
GenericType::GenericType(){
}

ostream& operator<<(ostream &stream, const GenericType& ref){
  ref->print(stream);
  return stream;
}

GenericType::GenericType(bool b){
  assignNode(new IntType(b));
}

GenericType::GenericType(int i){
  assignNode(new IntType(i));
}

GenericType::GenericType(double d){
  assignNode(new DoubleType(d));
}

GenericType::GenericType(const vector<int>& iv){
  assignNode(new IntVectorType(iv));
}

GenericType::GenericType(const vector<bool>& b_vec){
  vector<int> i_vec(b_vec.size());
  copy(b_vec.begin(),b_vec.end(), i_vec.begin());
  assignNode(new IntVectorType(i_vec));
}

GenericType::GenericType(const vector<double>& dv){
  assignNode(new DoubleVectorType(dv));
}

GenericType::GenericType(const string& s){
  assignNode(new StringType(s));
}

GenericType::GenericType(const char s[]){ 
  assignNode(new StringType(s));
}

bool GenericType::toBool() const{
  return bool(toInt());
}

int GenericType::toInt() const{
  if(isDouble()){
    double v = toDouble();
    casadi_assert_message(v == std::floor(v),"The value is not an integer");
    return int(v);
  } else {
    casadi_assert_message(isInt(),"type mismatch");
    return static_cast<const IntType*>(get())->d_;
  }
}

double GenericType::toDouble() const{
  if(isInt()){
    return double(toInt());
  } else {
    casadi_assert_message(isDouble(),"type mismatch");
    return static_cast<const DoubleType*>(get())->d_;
  }
}

const string& GenericType::toString() const{
  casadi_assert_message(isString(),"type mismatch");
  return static_cast<const StringType*>(get())->d_;
}

const vector<int>& GenericType::toIntVector() const{
  casadi_assert_message(isIntVector(),"type mismatch");
  return static_cast<const IntVectorType*>(get())->d_;
}

const vector<double>& GenericType::toDoubleVector() const{
  casadi_assert_message(isDoubleVector(),"type mismatch");
  return static_cast<const DoubleVectorType*>(get())->d_;
}


bool GenericType::operator==(const GenericType& op2) const{
  return !(*this != op2);
}

bool GenericType::operator!=(const GenericType& op2) const{
  if(isString() && op2.isString()){
    return toString().compare(op2.toString()) != 0;
  }

  if(isInt() && op2.isInt()){
    return toInt() != op2.toInt();
  }

  if(isDouble() && op2.isDouble()){
    return toDouble() != op2.toDouble();
  }

  if(isDoubleVector() && op2.isDoubleVector()){
    const vector<double> &v1 = toDoubleVector();
    const vector<double> &v2 = op2.toDoubleVector();
    if(v1.size() != v2.size()) return true;
    for(int i=0; i<v1.size(); ++i)
      if(v1[i] != v2[i]) return true;
    return false;
  }

  if(isIntVector() && op2.isIntVector()){
    const vector<int> &v1 = toIntVector();
    const vector<int> &v2 = op2.toIntVector();
    if(v1.size() != v2.size()) return true;
    for(int i=0; i<v1.size(); ++i)
      if(v1[i] != v2[i]) return true;
    return false;
  }
  
  // Different types
  return true;
}



} // namespace CasADi

