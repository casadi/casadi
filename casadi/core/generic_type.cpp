/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
#include "std_vector_tools.hpp"
#include "casadi_exception.hpp"
#include <cmath>

#include "function/function.hpp"
#include "functor.hpp"

using namespace std;

namespace casadi {

/// \cond INTERNAL
  typedef GenericTypeInternal<std::string> StringType;
  typedef GenericTypeInternal<double> DoubleType;
  typedef GenericTypeInternal<int> IntType;
  typedef GenericTypeInternal<bool> BoolType;
  typedef GenericTypeInternal<std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<std::vector<int> > IntVectorType;
  typedef GenericTypeInternal<std::vector<std::string> > StringVectorType;
  typedef GenericTypeInternal<SharedObject> SharedObjectType;
  typedef GenericTypeInternal<Function> FunctionType;
  typedef GenericTypeInternal<Dictionary> DictionaryType;
/// \endcond

opt_type GenericType::getType() const {
  return type_;
}

bool GenericType::can_cast_to(opt_type other) const {
  switch (other) {
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
  switch (type) {
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
  switch (type) {
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
      case OT_DERIVATIVEGENERATOR:
              return "OT_DERIVATIVEGENERATOR";
      case OT_CALLBACK:
              return "OT_CALLBACK";
      case OT_FUNCTION:
              return "OT_FUNCTION";
      case OT_VOIDPTR:
              return "OT_VOIDPTR";
      default:
              return "OT_UNKNOWN";

    }
}


bool GenericType::isBool() const {
  return is_a<bool>();
}

bool GenericType::isInt() const {
  return is_a<int>();
}

bool GenericType::isDouble() const {
  return is_a<double>();
}

bool GenericType::isString() const {
  return is_a<string>();
}

bool GenericType::isEmptyVector() const {
  return (isIntVector() && toIntVector().size()==0) ||
      (isDoubleVector() && toDoubleVector().size()==0) ||
      (isStringVector() && toStringVector().size()==0);
}

bool GenericType::isIntVector() const {
  return is_a<vector<int> >();
}

bool GenericType::isDoubleVector() const {
  return is_a<vector<double> >();
}

bool GenericType::isStringVector() const {
  return is_a<vector<string> >();
}

bool GenericType::isSharedObject() const {
  return is_a<SharedObject>();
}

bool GenericType::isFunction() const {
  return is_a<Function>();
}

bool GenericType::isDictionary() const {
  return is_a<Dictionary>();
}


GenericType::GenericType() {
}


ostream& operator<<(ostream &stream, const GenericType& ref) {
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

GenericType::GenericType(const vector<int>& iv) : type_(OT_INTEGERVECTOR) {
  assignNode(new IntVectorType(iv));
}

GenericType::GenericType(const vector<bool>& b_vec) : type_(OT_BOOLVECTOR) {
  vector<int> i_vec(b_vec.size());
  copy(b_vec.begin(), b_vec.end(), i_vec.begin());
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

GenericType::GenericType(const Function& f) : type_(OT_FUNCTION) {
  assignNode(new FunctionType(f));
}

  GenericType::GenericType(const DerivativeGenerator& f) : type_(OT_DERIVATIVEGENERATOR) {
    assignNode(new GenericTypeInternal<DerivativeGenerator>(f));
  }

GenericType::GenericType(const Callback& f) : type_(OT_CALLBACK) {
  assignNode(new GenericTypeInternal<Callback>(f));
}

bool GenericType::toBool() const {
  if (isBool()) {
    return static_cast<const BoolType*>(get())->d_;
  } else if (isInt()) {
    return static_cast<bool>(toInt());
  } else {
    casadi_assert_message(isBool(), "type mismatch");
    return false;
  }
}

int GenericType::toInt() const {
  if (isDouble()) {
    double v = toDouble();
    casadi_assert_message(v == std::floor(v), "The value is not an integer");
    return static_cast<int>(v);
  } else if (isBool()) {
    return static_cast<int>(toBool());
  } else {
    casadi_assert_message(isInt(), "type mismatch");
    return static_cast<const IntType*>(get())->d_;
  }
}

double GenericType::toDouble() const {
  if (isInt()) {
    return static_cast<double>(toInt());
  } else {
    casadi_assert_message(isDouble(), "type mismatch");
    return static_cast<const DoubleType*>(get())->d_;
  }
}

const string& GenericType::toString() const {
  casadi_assert_message(isString(), "type mismatch");
  return static_cast<const StringType*>(get())->d_;
}

const vector<int>& GenericType::toIntVector() const {
  casadi_assert_message(isIntVector(), "type mismatch");
  return static_cast<const IntVectorType*>(get())->d_;
}

const vector<double>& GenericType::toDoubleVector() const {
  casadi_assert_message(isDoubleVector(), "type mismatch");
  return static_cast<const DoubleVectorType*>(get())->d_;
}

vector<int>& GenericType::toIntVector() {
  casadi_assert_message(isIntVector(), "type mismatch");
  return static_cast<IntVectorType*>(get())->d_;
}

vector<double>& GenericType::toDoubleVector() {
  casadi_assert_message(isDoubleVector(), "type mismatch");
  return static_cast<DoubleVectorType*>(get())->d_;
}

const vector<string>& GenericType::toStringVector() const {
  casadi_assert_message(isStringVector(), "type mismatch");
  return static_cast<const StringVectorType*>(get())->d_;
}

const SharedObject& GenericType::toSharedObject() const {
  casadi_assert_message(isSharedObject(), "type mismatch");
  return static_cast<const SharedObjectType*>(get())->d_;
}

const Dictionary& GenericType::toDictionary() const {
  casadi_assert_message(isDictionary(), "type mismatch");
  return static_cast<const DictionaryType*>(get())->d_;
}

Dictionary& GenericType::toDictionary() {
  casadi_assert_message(isDictionary(), "type mismatch");
  return static_cast<DictionaryType*>(get())->d_;
}

const Function& GenericType::toFunction() const {
  casadi_assert_message(isFunction(), "type mismatch");
  return static_cast<const FunctionType*>(get())->d_;
}

bool GenericType::operator==(const GenericType& op2) const {
  return !(*this != op2);
}

bool GenericType::operator!=(const GenericType& op2) const {
  if (isString() && op2.isString()) {
    return toString().compare(op2.toString()) != 0;
  }

  if (isInt() && op2.isInt()) {
    return toInt() != op2.toInt();
  }

  if (isDouble() && op2.isDouble()) {
    return toDouble() != op2.toDouble();
  }

  if (isDoubleVector() && op2.isDoubleVector()) {
    const vector<double> &v1 = toDoubleVector();
    const vector<double> &v2 = op2.toDoubleVector();
    if (v1.size() != v2.size()) return true;
    for (int i=0; i<v1.size(); ++i)
      if (v1[i] != v2[i]) return true;
    return false;
  }

  if (isIntVector() && op2.isIntVector()) {
    const vector<int> &v1 = toIntVector();
    const vector<int> &v2 = op2.toIntVector();
    if (v1.size() != v2.size()) return true;
    for (int i=0; i<v1.size(); ++i)
      if (v1[i] != v2[i]) return true;
    return false;
  }

  // Different types
  return true;
}

GenericType::operator const DerivativeGenerator &() const {
    casadi_assert_message(is_a<DerivativeGenerator>(), "type mismatch");
    return static_cast<const GenericTypeInternal<DerivativeGenerator>*>(get())->d_;
  }

GenericType::operator const Callback &() const {
  casadi_assert_message(is_a<Callback>(), "type mismatch");
  return static_cast<const GenericTypeInternal<Callback>*>(get())->d_;
}

GenericType::GenericType(const Dictionary& dict) : type_(OT_DICTIONARY) {
  assignNode(new GenericTypeInternal<Dictionary>(dict));
}

GenericType::operator const GenericType::Dictionary& () const {
  casadi_assert_message(is_a<Dictionary>(), "type mismatch");
  return static_cast<const GenericTypeInternal<Dictionary>*>(get())->d_;
}

GenericType::operator GenericType::Dictionary& () {
  casadi_assert_message(is_a<Dictionary>(), "type mismatch");
  return static_cast<GenericTypeInternal<Dictionary>*>(get())->d_;
}

//GenericType::operator void*() const {
//  casadi_assert_message(is_a<void*>(), "type mismatch");
//  return static_cast<const GenericTypeInternal<void*>*>(get())->d_;
//}

void * GenericType::toVoidPointer() const {
  casadi_assert_message(is_a<void*>(), "type mismatch");
  return static_cast<const GenericTypeInternal<void*>*>(get())->d_;
}

GenericType::GenericType(void* ptr) : type_(OT_VOIDPTR) {
  assignNode(new GenericTypeInternal<void*>(ptr));
}


} // namespace casadi

