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

using namespace std;

namespace CasADi{


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

GenericTypeInternal* GenericType::operator->(){
  return (GenericTypeInternal*)SharedObject::operator->();
}

const GenericTypeInternal* GenericType::operator->() const{
  return (const GenericTypeInternal*)SharedObject::operator->();
}

bool GenericType::toBool() const{
  return bool((*this)->toInt());
}

int GenericType::toInt() const{
  return (*this)->toInt();  
}

double GenericType::toDouble() const{
  return (*this)->toDouble();    
}

const string& GenericType::toString() const{
  return (*this)->toString();    
}

const vector<int>& GenericType::toIntVector() const{
  return (*this)->toIntVector();    
}

const vector<double>& GenericType::toDoubleVector() const{
  return (*this)->toDoubleVector();      
}


bool GenericType::operator==(const GenericType& op2) const{
  return !(*this != op2);
}

bool GenericType::operator!=(const GenericType& op2) const{
  if(is_a<StringType>() && op2.is_a<StringType>()){
    return toString().compare(op2.toString()) != 0;
  }

  if(is_a<IntType>() && op2.is_a<IntType>()){
    return toInt() != op2.toInt();
  }

  if(is_a<DoubleType>() && op2.is_a<DoubleType>()){
    return toDouble() != op2.toDouble();
  }

  if(is_a<DoubleVectorType>() && op2.is_a<DoubleVectorType>()){
    const vector<double> &v1 = toDoubleVector();
    const vector<double> &v2 = op2.toDoubleVector();
    if(v1.size() != v2.size()) return true;
    for(int i=0; i<v1.size(); ++i)
      if(v1[i] != v2[i]) return true;
    return false;
  }

  if(is_a<IntVectorType>() && op2.is_a<IntVectorType>()){
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

