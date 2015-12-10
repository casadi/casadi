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

using namespace std;

namespace casadi {

  /// \cond INTERNAL
  typedef GenericTypeInternal<OT_STRING, std::string> StringType;
  typedef GenericTypeInternal<OT_REAL, double> DoubleType;
  typedef GenericTypeInternal<OT_INTEGER, int> IntType;
  typedef GenericTypeInternal<OT_BOOLEAN, bool> BoolType;
  typedef GenericTypeInternal<OT_REALVECTOR, std::vector<double> > DoubleVectorType;
  typedef GenericTypeInternal<OT_INTEGERVECTOR, std::vector<int> > IntVectorType;
  typedef GenericTypeInternal<OT_INTEGERVECTORVECTOR,
                              std::vector< std::vector<int> > > IntVectorVectorType;
  typedef GenericTypeInternal<OT_STRINGVECTOR, std::vector<std::string> > StringVectorType;
  typedef GenericTypeInternal<OT_FUNCTION, Function> FunctionType;
  typedef GenericTypeInternal<OT_DICT, Dict> DictType;
  typedef GenericTypeInternal<OT_VOIDPTR, void*> VoidPointerType;
  /// \endcond

  bool GenericType::can_cast_to(TypeID other) const {
    switch (other) {
    case OT_BOOLEAN:
      return isBool() || isInt() || isDouble();
    case OT_BOOLVECTOR:
      return isIntVector() || isDoubleVector();
    case OT_INTEGER:
    case OT_REAL:
      return isInt() || isDouble();
    case OT_INTEGERVECTOR:
    case OT_REALVECTOR:
      return isDoubleVector() || isIntVector();
    case OT_STRINGVECTOR:
      return isStringVector() || isString();
    default:
      return getType() == other;
    }
  }

  GenericType GenericType::from_type(TypeID type) {
    switch (type) {
    case OT_INTEGERVECTOR:
      return std::vector<int>();
    case OT_INTEGERVECTORVECTOR:
      return std::vector< std::vector<int> >();
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

  std::string GenericType::get_type_description(TypeID type) {
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
    case OT_INTEGERVECTORVECTOR:
      return "OT_INTEGERVECTORVECTOR";
    case OT_BOOLVECTOR:
      return "OT_BOOLVECTOR";
    case OT_REALVECTOR:
      return "OT_REALVECTOR";
    case OT_STRINGVECTOR:
      return "OT_STRINGVECTOR";
    case OT_DICT:
      return "OT_DICT";
    case OT_FUNCTION:
      return "OT_FUNCTION";
    case OT_VOIDPTR:
      return "OT_VOIDPTR";
    default:
      return "OT_UNKNOWN";

    }
  }


  bool GenericType::isBool() const {
    return getType()==OT_BOOLEAN;
  }

  bool GenericType::isInt() const {
    return getType()==OT_INTEGER;
  }

  bool GenericType::isDouble() const {
    return getType()==OT_REAL;
  }

  bool GenericType::isString() const {
    return getType()==OT_STRING;
  }

  bool GenericType::is_emptyVector() const {
    return (isIntVector() && toIntVector().size()==0) ||
      (isIntVectorVector() && toIntVectorVector().size()==0) ||
      (isDoubleVector() && toDoubleVector().size()==0) ||
      (isStringVector() && toStringVector().size()==0);
  }

  bool GenericType::isIntVector() const {
    return getType()==OT_INTEGERVECTOR;
  }

  bool GenericType::isIntVectorVector() const {
    return getType()==OT_INTEGERVECTORVECTOR;
  }

  bool GenericType::isDoubleVector() const {
    return getType()==OT_REALVECTOR;
  }

  bool GenericType::isStringVector() const {
    return getType()==OT_STRINGVECTOR;
  }

  bool GenericType::isFunction() const {
    return getType()==OT_FUNCTION;
  }

  bool GenericType::isVoidPointer() const {
    return getType()==OT_VOIDPTR;
  }

  bool GenericType::isDict() const {
    return getType()==OT_DICT;
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

  GenericType::GenericType(bool b) {
    assignNode(new BoolType(b));
  }

  GenericType::GenericType(int i) {
    assignNode(new IntType(i));
  }

  GenericType::GenericType(double d) {
    assignNode(new DoubleType(d));
  }

  GenericType::GenericType(const vector<int>& iv) {
    assignNode(new IntVectorType(iv));
  }

  GenericType::GenericType(const vector<vector<int> >& ivv) {
    assignNode(new IntVectorVectorType(ivv));
  }

  GenericType::GenericType(const vector<bool>& b_vec) {
    vector<int> i_vec(b_vec.size());
    copy(b_vec.begin(), b_vec.end(), i_vec.begin());
    assignNode(new IntVectorType(i_vec));
  }

  GenericType::GenericType(const vector<double>& dv) {
    assignNode(new DoubleVectorType(dv));
  }

  GenericType::GenericType(const vector<string>& sv) {
    assignNode(new StringVectorType(sv));
  }

  GenericType::GenericType(const string& s) {
    assignNode(new StringType(s));
  }

  GenericType::GenericType(const char s[]) {
    assignNode(new StringType(s));
  }

  GenericType::GenericType(const Function& f) {
    assignNode(new FunctionType(f));
  }

  const bool& GenericType::asBool() const {
    casadi_assert(isBool());
    return static_cast<const BoolType*>(get())->d_;
  }

  const int& GenericType::asInt() const {
    casadi_assert(isInt());
    return static_cast<const IntType*>(get())->d_;
  }

  const double& GenericType::asDouble() const {
    casadi_assert(isDouble());
    return static_cast<const DoubleType*>(get())->d_;
  }

  const std::string& GenericType::asString() const {
    casadi_assert(isString());
    return static_cast<const StringType*>(get())->d_;
  }

  const std::vector<int>& GenericType::asIntVector() const {
    casadi_assert(isIntVector());
    return static_cast<const IntVectorType*>(get())->d_;
  }

  const std::vector<std::vector<int> >& GenericType::asIntVectorVector() const {
    casadi_assert(isIntVectorVector());
    return static_cast<const IntVectorVectorType*>(get())->d_;
  }

  const std::vector<double>& GenericType::asDoubleVector() const {
    casadi_assert(isDoubleVector());
    return static_cast<const DoubleVectorType*>(get())->d_;
  }

  const std::vector<std::string>& GenericType::asStringVector() const {
    casadi_assert(isStringVector());
    return static_cast<const StringVectorType*>(get())->d_;
  }

  const GenericType::Dict& GenericType::asDict() const {
    casadi_assert(isDict());
    return static_cast<const DictType*>(get())->d_;
  }

  const Function& GenericType::asFunction() const {
    casadi_assert(isFunction());
    return static_cast<const FunctionType*>(get())->d_;
  }

  void* const & GenericType::asVoidPointer() const {
    casadi_assert(isVoidPointer());
    return static_cast<const VoidPointerType*>(get())->d_;
  }

  bool GenericType::toBool() const {
    if (isBool()) {
      return asBool();
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
      return asInt();
    }
  }

  double GenericType::toDouble() const {
    if (isInt()) {
      return static_cast<double>(toInt());
    } else {
      casadi_assert_message(isDouble(), "type mismatch");
      return asDouble();
    }
  }

  string GenericType::toString() const {
    casadi_assert_message(isString(), "type mismatch");
    return asString();
  }

  vector<int> GenericType::toIntVector() const {
    casadi_assert_message(isIntVector(), "type mismatch");
    return asIntVector();
  }

  vector<vector<int> > GenericType::toIntVectorVector() const {
    casadi_assert_message(isIntVectorVector(), "type mismatch");
    return asIntVectorVector();
  }

  vector<double> GenericType::toDoubleVector() const {
    casadi_assert_message(isDoubleVector(), "type mismatch");
    return asDoubleVector();
  }

  vector<string> GenericType::toStringVector() const {
    if (isString()) {
      std::string s = asString();
      return vector<string>(1, s);
    } else {
      casadi_assert_message(isStringVector(), "type mismatch");
      return asStringVector();
    }
  }

  Dict GenericType::toDict() const {
    casadi_assert_message(isDict(), "type mismatch");
    return asDict();
  }

  Function GenericType::toFunction() const {
    casadi_assert_message(isFunction(), "type mismatch");
    return asFunction();
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

    if (isIntVectorVector() && op2.isIntVectorVector()) {
      const vector< vector<int> > &v1 = toIntVectorVector();
      const vector< vector<int> > &v2 = op2.toIntVectorVector();
      if (v1.size() != v2.size()) return true;
      for (int i=0; i<v1.size(); ++i) {
        if (v1[i].size() != v2[i].size()) return true;
        for (int j=0; j<v1[i].size(); ++j) {
          if (v1[i][j] != v2[i][j]) return true;
        }
      }
      return false;
    }

    // Different types
    return true;
  }

  GenericType::GenericType(const Dict& dict) {
    assignNode(new DictType(dict));
  }

  void* GenericType::toVoidPointer() const {
    casadi_assert_message(getType()==OT_VOIDPTR, "type mismatch");
    return asVoidPointer();
  }

  GenericType::GenericType(void* ptr) {
    assignNode(new VoidPointerType(ptr));
  }

  TypeID GenericType::getType() const {
    if (isNull()) {
      return OT_NULL;
    } else {
      return static_cast<const GenericTypeBase*>(get())->getType();
    }
  }

} // namespace casadi
