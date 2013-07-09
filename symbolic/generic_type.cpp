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

#include "fx/fx.hpp"

using namespace std;

namespace CasADi{

  typedef GenericTypeInternal<std::string> StringType;
  typedef GenericTypeInternal<double> DoubleType;
  typedef GenericTypeInternal<int> IntType;
  typedef GenericTypeInternal<bool> BoolType;
  typedef GenericTypeInternal<std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<std::vector<int> > IntVectorType;
  typedef GenericTypeInternal<std::vector<std::string> > StringVectorType;
  typedef GenericTypeInternal<SharedObject> SharedObjectType;
  typedef GenericTypeInternal<FX> FXType;
  typedef GenericTypeInternal<Dictionary> DictionaryType;
  
opt_type GenericType::getType() const {
  return type_;
}
    
bool GenericType::can_cast_to(opt_type other) const {
  switch(other)
    {
      case OT_BOOLEAN:
        return isBool() || isInt() || isDouble();
      case OT_BOOLVECTOR:
        return isIntVector() || isDoubleVector();
      case OT_INTEGER: case OT_REAL:
        return isInt() || isDouble();
      case OT_INTEGERVECTOR: case OT_REALVECTOR:
        return isDoubleVector() || isIntVector();
      default:
        return type_ == other;
  }
}

GenericType GenericType::from_type(opt_type type) {
  switch(type)
    {
      case OT_INTEGERVECTOR:
              return std::vector<int>();
      case OT_BOOLVECTOR:
              return std::vector<bool>();
      case OT_REALVECTOR:
              return std::vector<double>();
      case OT_STRINGVECTOR:
              return std::vector<std::string>();
      default:
              casadi_error("empty_from_type. Unsupported type " << type);
    }
}

std::string GenericType::get_type_description(const opt_type &type) {
  switch(type)
    {
      case OT_BOOLEAN:
              return "OT_BOOLEAN";
      case OT_INTEGER:
              return "OT_INTEGER";
      case OT_REAL:
              return "OT_REAL";
      case OT_STRING:
              return "OT_STRING";
      case OT_INTEGERVECTOR:
              return "OT_INTEGERVECTOR";
      case OT_BOOLVECTOR:
              return "OT_BOOLVECTOR";
      case OT_REALVECTOR:
              return "OT_REALVECTOR";
      case OT_STRINGVECTOR:
              return "OT_STRINGVECTOR";
      case OT_DICTIONARY:
              return "OT_DICTIONARY";
      case OT_NLPSOLVER:
              return "OT_NLPSOLVER";
      case OT_LINEARSOLVER:
              return "OT_LINEARSOLVER";
      case OT_LPSOLVER:
              return "OT_LPSOLVER";
      case OT_INTEGRATOR:
              return "OT_INTEGRATOR";
      case OT_QPSOLVER:
              return "OT_QPSOLVER";
      case OT_QCQPSOLVER:
              return "OT_QCQPSOLVER";
      case OT_SOCPSOLVER:
              return "OT_SOCPSOLVER";
      case OT_SDPSOLVER:
              return "OT_SDPSOLVER";
      case OT_SDQPSOLVER:
              return "OT_SDQPSOLVER";
      case OT_IMPLICITFUNCTION:
              return "OT_IMPLICITFUNCTION";
      case OT_JACOBIANGENERATOR:
              return "OT_JACOBIANGENERATOR";
      case OT_SPARSITYGENERATOR:
              return "OT_SPARSITYGENERATOR";
      case OT_FX:
              return "OT_FX";
      case OT_VOIDPTR:
              return "OT_VOIDPTR";
      default:
              return "OT_UNKNOWN";
              
    }
};
    

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

bool GenericType::isEmptyVector() const{
  return (isIntVector() && toIntVector().size()==0 ) || (isDoubleVector() && toDoubleVector().size()==0 ) || (isStringVector() && toStringVector().size()==0 );
}

bool GenericType::isIntVector() const{
  return is_a<vector<int> >();
}
    
bool GenericType::isDoubleVector() const{
  return is_a<vector<double> >();
}
  
bool GenericType::isStringVector() const{
  return is_a<vector<string> >();
}
  
bool GenericType::isSharedObject() const{
  return is_a<SharedObject>();
}
  
bool GenericType::isFX() const{
  return is_a<FX>();
}

bool GenericType::isDictionary() const{
  return is_a<Dictionary>();
}


GenericType::GenericType(){
}


ostream& operator<<(ostream &stream, const GenericType& ref){
  if (ref.isNull()) {
    stream << "None";
  } else {
    ref->print(stream);
  }
  return stream;
}

GenericType::GenericType(bool b) : type_(OT_BOOLEAN) {
  assignNode(new BoolType(b));
}

GenericType::GenericType(int i) : type_(OT_INTEGER) {
  assignNode(new IntType(i));
}

GenericType::GenericType(double d) : type_(OT_REAL) {
  assignNode(new DoubleType(d));
}

GenericType::GenericType(const vector<int>& iv) : type_(OT_INTEGERVECTOR){
  assignNode(new IntVectorType(iv));
}

GenericType::GenericType(const vector<bool>& b_vec) : type_(OT_BOOLVECTOR){
  vector<int> i_vec(b_vec.size());
  copy(b_vec.begin(),b_vec.end(), i_vec.begin());
  assignNode(new IntVectorType(i_vec));
}

GenericType::GenericType(const vector<double>& dv) : type_(OT_REALVECTOR) {
  assignNode(new DoubleVectorType(dv));
}

GenericType::GenericType(const vector<string>& sv) : type_(OT_STRINGVECTOR) {
  assignNode(new StringVectorType(sv));
}

GenericType::GenericType(const string& s) : type_(OT_STRING) {
  assignNode(new StringType(s));
}

GenericType::GenericType(const char s[])  : type_(OT_STRING) { 
  assignNode(new StringType(s));
}

//GenericType::GenericType(const GenericType& obj) {
//  type_ = obj.type_;
//  assignNode(new SharedObjectType(obj));
//}

GenericType::GenericType(const SharedObject& obj) : type_(OT_UNKNOWN) {
  assignNode(new SharedObjectType(obj));
}

GenericType::GenericType(const FX& f) : type_(OT_FX) {
  assignNode(new FXType(f));
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

const vector<string>& GenericType::toStringVector() const{
  casadi_assert_message(isStringVector(),"type mismatch");
  return static_cast<const StringVectorType*>(get())->d_;
}

const SharedObject& GenericType::toSharedObject() const{
  casadi_assert_message(isSharedObject(),"type mismatch");
  return static_cast<const SharedObjectType*>(get())->d_;
}

const Dictionary& GenericType::toDictionary() const{
  casadi_assert_message(isDictionary(),"type mismatch");
  return static_cast<const DictionaryType*>(get())->d_;
}

const FX& GenericType::toFX() const{
  casadi_assert_message(isFX(),"type mismatch");
  return static_cast<const FXType*>(get())->d_;
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
  
GenericType::GenericType(NLPSolverCreator ptr)  : type_(OT_NLPSOLVER) {
  assignNode(new GenericTypeInternal<NLPSolverCreator>(ptr));
}

GenericType::GenericType(linearSolverCreator ptr) : type_(OT_LINEARSOLVER) {
  assignNode(new GenericTypeInternal<linearSolverCreator >(ptr));
}

GenericType::GenericType(integratorCreator ptr) : type_(OT_INTEGRATOR) {
  assignNode(new GenericTypeInternal<integratorCreator>(ptr));
}

GenericType::GenericType(QPSolverCreator ptr) : type_(OT_QPSOLVER) {
  assignNode(new GenericTypeInternal<QPSolverCreator>(ptr));
}

GenericType::GenericType(LPSolverCreator ptr) : type_(OT_LPSOLVER) {
  assignNode(new GenericTypeInternal<LPSolverCreator>(ptr));
}

GenericType::GenericType(SDPSolverCreator ptr) : type_(OT_SDPSOLVER) {
  assignNode(new GenericTypeInternal<SDPSolverCreator>(ptr));
}

GenericType::GenericType(SDQPSolverCreator ptr) : type_(OT_SDQPSOLVER) {
  assignNode(new GenericTypeInternal<SDQPSolverCreator>(ptr));
}

GenericType::GenericType(SOCPSolverCreator ptr) : type_(OT_SOCPSOLVER) {
  assignNode(new GenericTypeInternal<SOCPSolverCreator>(ptr));
}

GenericType::GenericType(QCQPSolverCreator ptr) : type_(OT_QCQPSOLVER) {
  assignNode(new GenericTypeInternal<QCQPSolverCreator>(ptr));
}

GenericType::GenericType(implicitFunctionCreator ptr) : type_(OT_IMPLICITFUNCTION) {
  assignNode(new GenericTypeInternal<implicitFunctionCreator>(ptr));
}

GenericType::GenericType(JacobianGenerator ptr) : type_(OT_JACOBIANGENERATOR) {
  assignNode(new GenericTypeInternal<JacobianGenerator>(ptr));
}

GenericType::GenericType(SparsityGenerator ptr) : type_(OT_SPARSITYGENERATOR) {
  assignNode(new GenericTypeInternal<SparsityGenerator>(ptr));
}

GenericType::operator NLPSolverCreator() const{
  casadi_assert_message(is_a<NLPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<NLPSolverCreator>*>(get())->d_;
}

GenericType::operator LPSolverCreator() const{
  casadi_assert_message(is_a<LPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<LPSolverCreator>*>(get())->d_;
}

GenericType::operator SDPSolverCreator() const{
  casadi_assert_message(is_a<SDPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<SDPSolverCreator>*>(get())->d_;
}

GenericType::operator SDQPSolverCreator() const{
  casadi_assert_message(is_a<SDQPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<SDQPSolverCreator>*>(get())->d_;
}

GenericType::operator SOCPSolverCreator() const{
  casadi_assert_message(is_a<SOCPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<SOCPSolverCreator>*>(get())->d_;
}

GenericType::operator QCQPSolverCreator() const{
  casadi_assert_message(is_a<QCQPSolverCreator>(),"type mismatch");
  return static_cast<const GenericTypeInternal<QCQPSolverCreator>*>(get())->d_;
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

GenericType::GenericType(const Dictionary& dict) : type_(OT_DICTIONARY) {
  assignNode(new GenericTypeInternal<Dictionary>(dict));
}

GenericType::operator const GenericType::Dictionary& () const{
  casadi_assert_message(is_a<Dictionary>(),"type mismatch");
  return static_cast<const GenericTypeInternal<Dictionary>*>(get())->d_;
}

//GenericType::operator void*() const{
//  casadi_assert_message(is_a<void*>(),"type mismatch");
//  return static_cast<const GenericTypeInternal<void*>*>(get())->d_;
//}

void * GenericType::toVoidPointer() const {
  casadi_assert_message(is_a<void*>(),"type mismatch");
  return static_cast<const GenericTypeInternal<void*>*>(get())->d_; 
}

GenericType::GenericType(void* ptr) : type_(OT_VOIDPTR){
  assignNode(new GenericTypeInternal<void*>(ptr));
}


} // namespace CasADi

