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
  typedef GenericTypeInternal<bool> BoolType;
  typedef GenericTypeInternal<std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<std::vector<int> > IntVectorType;
  typedef GenericTypeInternal<SharedObject> SharedObjectType;

bool GenericType::isBool() const{
  return is_a<bool>();
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
  
bool GenericType::isSharedObject() const{
  return is_a<SharedObject>();
}
  
GenericType::GenericType(){
}

ostream& operator<<(ostream &stream, const GenericType& ref){
  ref->print(stream);
  return stream;
}

GenericType::GenericType(bool b){
  assignNode(new BoolType(b));
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

GenericType::GenericType(const SharedObject& obj){
  assignNode(new SharedObjectType(obj));
}


bool GenericType::toBool() const{
  if (isBool()) {
    return static_cast<const BoolType*>(get())->d_;
  } else if (isInt()) {
    return bool(toInt());
  } else {
	casadi_assert_message(isBool(),"type mismatch");
	return false;
  }
}

int GenericType::toInt() const{
  if(isDouble()){
    double v = toDouble();
    casadi_assert_message(v == std::floor(v),"The value is not an integer");
    return int(v);
  } else if (isBool()) {
    return int(toBool());
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

const SharedObject& GenericType::toSharedObject() const{
  casadi_assert_message(isSharedObject(),"type mismatch");
  return static_cast<const SharedObjectType*>(get())->d_;
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

GenericType::GenericType(NLPSolverCreator ptr){
  assignNode(new GenericTypeInternal<NLPSolverCreator>(ptr));
}

GenericType::GenericType(linearSolverCreator ptr){
  assignNode(new GenericTypeInternal<linearSolverCreator >(ptr));
}

GenericType::GenericType(integratorCreator ptr){
  assignNode(new GenericTypeInternal<integratorCreator>(ptr));
}

GenericType::GenericType(QPSolverCreator ptr){
  assignNode(new GenericTypeInternal<QPSolverCreator>(ptr));
}

GenericType::GenericType(implicitFunctionCreator ptr){
  assignNode(new GenericTypeInternal<implicitFunctionCreator>(ptr));
}

GenericType::GenericType(JacobianGenerator ptr){
  assignNode(new GenericTypeInternal<JacobianGenerator>(ptr));
}

GenericType::GenericType(SparsityGenerator ptr){
  assignNode(new GenericTypeInternal<SparsityGenerator>(ptr));
}

GenericType::operator NLPSolverCreator() const{
  casadi_assert_message(is_a<NLPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<NLPSolverCreator>*>(get())->d_;
}

GenericType::operator linearSolverCreator() const{
  casadi_assert_message(is_a<linearSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<linearSolverCreator>*>(get())->d_;
}
GenericType::operator integratorCreator() const{
  casadi_assert_message(is_a<integratorCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<integratorCreator>*>(get())->d_;
}

GenericType::operator QPSolverCreator() const{
  casadi_assert_message(is_a<QPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<QPSolverCreator>*>(get())->d_;
}

GenericType::operator implicitFunctionCreator() const{
  casadi_assert_message(is_a<implicitFunctionCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<implicitFunctionCreator>*>(get())->d_;
}

GenericType::operator JacobianGenerator() const{
  casadi_assert_message(is_a<JacobianGenerator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<JacobianGenerator>*>(get())->d_;
}

GenericType::operator SparsityGenerator() const{
  casadi_assert_message(is_a<SparsityGenerator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<SparsityGenerator>*>(get())->d_;
}

GenericType::GenericType(const Dictionary& dict){
  assignNode(new GenericTypeInternal<Dictionary>(dict));
}

GenericType::operator const GenericType::Dictionary& () const{
  casadi_assert_message(is_a<Dictionary>(),"type mismatch");
  return static_cast<const GenericTypeInternal<Dictionary>*>(get())->d_;
}

GenericType::operator void*() const{
  casadi_assert_message(is_a<void*>(),"type mismatch");
  return static_cast<const GenericTypeInternal<void*>*>(get())->d_;
}

GenericType::GenericType(void* ptr){
  assignNode(new GenericTypeInternal<void*>(ptr));
}


} // namespace CasADi

