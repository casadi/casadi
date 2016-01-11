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
#include "exception.hpp"
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
      return is_bool() || is_int() || is_double();
    case OT_BOOLVECTOR:
      return is_int_vector() || is_double_vector();
    case OT_INTEGER:
    case OT_REAL:
      return is_int() || is_double();
    case OT_INTEGERVECTOR:
    case OT_REALVECTOR:
      return is_double_vector() || is_int_vector();
    case OT_STRINGVECTOR:
      return is_string_vector() || is_string();
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


  bool GenericType::is_bool() const {
    return getType()==OT_BOOLEAN;
  }

  bool GenericType::is_int() const {
    return getType()==OT_INTEGER;
  }

  bool GenericType::is_double() const {
    return getType()==OT_REAL;
  }

  bool GenericType::is_string() const {
    return getType()==OT_STRING;
  }

  bool GenericType::is_empty_vector() const {
    return (is_int_vector() && toIntVector().size()==0) ||
      (is_int_vector_vector() && toIntVectorVector().size()==0) ||
      (is_double_vector() && toDoubleVector().size()==0) ||
      (is_string_vector() && toStringVector().size()==0);
  }

  bool GenericType::is_int_vector() const {
    return getType()==OT_INTEGERVECTOR;
  }

  bool GenericType::is_int_vector_vector() const {
    return getType()==OT_INTEGERVECTORVECTOR;
  }

  bool GenericType::is_double_vector() const {
    return getType()==OT_REALVECTOR;
  }

  bool GenericType::is_string_vector() const {
    return getType()==OT_STRINGVECTOR;
  }

  bool GenericType::is_function() const {
    return getType()==OT_FUNCTION;
  }

  bool GenericType::is_void_pointer() const {
    return getType()==OT_VOIDPTR;
  }

  bool GenericType::is_dict() const {
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
    casadi_assert(is_bool());
    return static_cast<const BoolType*>(get())->d_;
  }

  const int& GenericType::asInt() const {
    casadi_assert(is_int());
    return static_cast<const IntType*>(get())->d_;
  }

  const double& GenericType::asDouble() const {
    casadi_assert(is_double());
    return static_cast<const DoubleType*>(get())->d_;
  }

  const std::string& GenericType::asString() const {
    casadi_assert(is_string());
    return static_cast<const StringType*>(get())->d_;
  }

  const std::vector<int>& GenericType::asIntVector() const {
    casadi_assert(is_int_vector());
    return static_cast<const IntVectorType*>(get())->d_;
  }

  const std::vector<std::vector<int> >& GenericType::asIntVectorVector() const {
    casadi_assert(is_int_vector_vector());
    return static_cast<const IntVectorVectorType*>(get())->d_;
  }

  const std::vector<double>& GenericType::asDoubleVector() const {
    casadi_assert(is_double_vector());
    return static_cast<const DoubleVectorType*>(get())->d_;
  }

  const std::vector<std::string>& GenericType::asStringVector() const {
    casadi_assert(is_string_vector());
    return static_cast<const StringVectorType*>(get())->d_;
  }

  const GenericType::Dict& GenericType::asDict() const {
    casadi_assert(is_dict());
    return static_cast<const DictType*>(get())->d_;
  }

  const Function& GenericType::asFunction() const {
    casadi_assert(is_function());
    return static_cast<const FunctionType*>(get())->d_;
  }

  void* const & GenericType::asVoidPointer() const {
    casadi_assert(is_void_pointer());
    return static_cast<const VoidPointerType*>(get())->d_;
  }

  bool GenericType::toBool() const {
    if (is_bool()) {
      return asBool();
    } else if (is_int()) {
      return static_cast<bool>(toInt());
    } else {
      casadi_assert_message(is_bool(), "type mismatch");
      return false;
    }
  }

  int GenericType::toInt() const {
    if (is_double()) {
      double v = toDouble();
      casadi_assert_message(v == std::floor(v), "The value is not an integer");
      return static_cast<int>(v);
    } else if (is_bool()) {
      return static_cast<int>(toBool());
    } else {
      casadi_assert_message(is_int(), "type mismatch");
      return asInt();
    }
  }

  double GenericType::toDouble() const {
    if (is_int()) {
      return static_cast<double>(toInt());
    } else {
      casadi_assert_message(is_double(), "type mismatch");
      return asDouble();
    }
  }

  string GenericType::toString() const {
    casadi_assert_message(is_string(), "type mismatch");
    return asString();
  }

  vector<int> GenericType::toIntVector() const {
    casadi_assert_message(is_int_vector(), "type mismatch");
    return asIntVector();
  }

  vector<vector<int> > GenericType::toIntVectorVector() const {
    casadi_assert_message(is_int_vector_vector(), "type mismatch");
    return asIntVectorVector();
  }

  vector<double> GenericType::toDoubleVector() const {
    casadi_assert_message(is_double_vector(), "type mismatch");
    return asDoubleVector();
  }

  vector<string> GenericType::toStringVector() const {
    if (is_string()) {
      std::string s = asString();
      return vector<string>(1, s);
    } else {
      casadi_assert_message(is_string_vector(), "type mismatch");
      return asStringVector();
    }
  }

  Dict GenericType::toDict() const {
    casadi_assert_message(is_dict(), "type mismatch");
    return asDict();
  }

  Function GenericType::toFunction() const {
    casadi_assert_message(is_function(), "type mismatch");
    return asFunction();
  }

  bool GenericType::operator==(const GenericType& op2) const {
    return !(*this != op2);
  }

  bool GenericType::operator!=(const GenericType& op2) const {
    if (is_string() && op2.is_string()) {
      return toString().compare(op2.toString()) != 0;
    }

    if (is_int() && op2.is_int()) {
      return toInt() != op2.toInt();
    }

    if (is_double() && op2.is_double()) {
      return toDouble() != op2.toDouble();
    }

    if (is_double_vector() && op2.is_double_vector()) {
      const vector<double> &v1 = toDoubleVector();
      const vector<double> &v2 = op2.toDoubleVector();
      if (v1.size() != v2.size()) return true;
      for (int i=0; i<v1.size(); ++i)
        if (v1[i] != v2[i]) return true;
      return false;
    }

    if (is_int_vector() && op2.is_int_vector()) {
      const vector<int> &v1 = toIntVector();
      const vector<int> &v2 = op2.toIntVector();
      if (v1.size() != v2.size()) return true;
      for (int i=0; i<v1.size(); ++i)
        if (v1[i] != v2[i]) return true;
      return false;
    }

    if (is_int_vector_vector() && op2.is_int_vector_vector()) {
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
