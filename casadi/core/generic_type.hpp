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


#ifndef CASADI_GENERIC_TYPE_HPP
#define CASADI_GENERIC_TYPE_HPP

#include "shared_object.hpp"
#include "casadi_types.hpp"
#include <string>
#include <vector>

namespace casadi {

  /** \brief  Types of options */
  enum TypeID {
    OT_NULL,
    OT_BOOLEAN,
    OT_INTEGER,
    OT_REAL,
    OT_STRING,
    OT_INTEGERVECTOR,
    OT_INTEGERVECTORVECTOR,
    OT_BOOLVECTOR,
    OT_REALVECTOR,
    OT_STRINGVECTOR,
    OT_DICT,
    OT_DERIVATIVEGENERATOR,
    OT_FUNCTION,
    OT_CALLBACK,
    OT_VOIDPTR,
    OT_UNKNOWN};

  /** \brief Generic data type, can hold different types such as bool, int, string etc.
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT GenericType : public SharedObject {
  public:
    /// C++ equivalent of Python's dict or MATLAB's struct
    typedef std::map<std::string, GenericType> Dict;

    /// Default constructor
    GenericType();

    /// Constructors (implicit type conversion)
    GenericType(bool b);
    GenericType(int i);
    GenericType(double d);
    GenericType(const std::string& s);
    GenericType(const std::vector<bool>& iv);
    GenericType(const std::vector<int>& iv);
    GenericType(const std::vector< std::vector<int> >& ivv);
    GenericType(const std::vector<double>& dv);
    GenericType(const std::vector<std::string>& sv);
    GenericType(const char s[]);
    GenericType(const Function& f);
    GenericType(const DerivativeGenerator& c);
    GenericType(const Callback& c);
    GenericType(const Dict& dict);
    #ifndef SWIG
    GenericType(void* ptr);
    #endif // SWIG

    /// Get a description of a type
    static std::string get_type_description(TypeID type);

    /// Get a description of the object's type
    std::string get_description() const { return get_type_description(getType()); }

    /// Construct a GenericType given an TypeID
    static GenericType from_type(TypeID type);

    #ifndef SWIG
    ///@{
    /// Implicit typecasting
    operator bool() const { return toBool();}
    operator int() const { return toInt();}
    operator double() const { return toDouble();}
    operator std::string() const { return toString();}
    operator std::vector<int>() const { return toIntVector();}
    operator std::vector<std::vector<int> >() const { return toIntVectorVector();}
    operator std::vector<double>() const { return toDoubleVector();}
    operator std::vector<std::string>() const { return toStringVector();}
    operator const Function&() const { return asFunction();}
    operator const Dict&() const { return asDict();}
    operator const DerivativeGenerator& () const { return asDerivativeGenerator();}
    operator const Callback& () const { return asCallback();}
    ///@}
    #endif // SWIG

    // Get type of object
    TypeID getType() const;

    bool can_cast_to(TypeID other) const;
    bool can_cast_to(const GenericType& other) const { return can_cast_to(other.getType()) ;}

    ///@{
    /** \brief Check if a particular type */
    bool isBool() const;
    bool isInt() const;
    bool isDouble() const;
    bool isString() const;
    bool isemptyVector() const;
    bool isIntVector() const;
    bool isIntVectorVector() const;
    bool isDoubleVector() const;
    bool isStringVector() const;
    bool isDict() const;
    bool isFunction() const;
    bool isVoidPointer() const;
    bool isCallback() const;
    bool isDerivativeGenerator() const;
    ///@}

#ifndef SWIG
    ///@{
    /** \brief Cast to the internal type */
    const bool& asBool() const;
    const int& asInt() const;
    const double& asDouble() const;
    const std::string& asString() const;
    const std::vector<int>& asIntVector() const;
    const std::vector<std::vector<int> >& asIntVectorVector() const;
    const std::vector<double>& asDoubleVector() const;
    const std::vector<std::string>& asStringVector() const;
    const Dict& asDict() const;
    const Function& asFunction() const;
    void* const & asVoidPointer() const;
    const DerivativeGenerator& asDerivativeGenerator() const;
    const Callback& asCallback() const;
    ///@}
#endif // SWIG

    ///@{
    //! \brief Convert to a type
    bool toBool() const;
    int toInt() const;
    double toDouble() const;
    std::string toString() const;
    std::vector<int> toIntVector() const;
    std::vector< std::vector<int> > toIntVectorVector() const;
    std::vector<double> toDoubleVector() const;
    std::vector<std::string> toStringVector() const;
    Dict toDict() const;
    Function toFunction() const;
    void* toVoidPointer() const;
    ///@}

    //! \brief Equality
    bool operator==(const GenericType& op2) const;
    bool operator!=(const GenericType& op2) const;

#ifndef SWIG
    //! \brief Print
    CASADI_EXPORT friend std::ostream& operator<<(std::ostream &stream,
                                                  const GenericType& ref);
#endif // SWIG
  };

  /// C++ equivalent of Python's dict or MATLAB's struct
  typedef GenericType::Dict Dict;

#ifndef SWIG
  // Create dictionary with 1 element
  inline Dict
  make_dict(const std::string& n0, const GenericType& x0) {
    Dict ret;
    ret[n0]=x0;
    return ret;
  }

  // Create dictionary with 2 elements
  inline Dict make_dict(const std::string& n0, const GenericType& x0,
                        const std::string& n1, const GenericType& x1) {
    Dict ret=make_dict(n0, x0);
    ret[n1]=x1;
    return ret;
  }

  // Create dictionary with 3 elements
  inline Dict make_dict(const std::string& n0, const GenericType& x0,
                        const std::string& n1, const GenericType& x1,
                        const std::string& n2, const GenericType& x2) {
    Dict ret=make_dict(n0, x0, n1, x1);
    ret[n2]=x2;
    return ret;
  }

  // Create dictionary with 4 elements
  inline Dict make_dict(const std::string& n0, const GenericType& x0,
                        const std::string& n1, const GenericType& x1,
                        const std::string& n2, const GenericType& x2,
                        const std::string& n3, const GenericType& x3) {
    Dict ret=make_dict(n0, x0, n1, x1, n2, x2);
    ret[n3]=x3;
    return ret;
  }

  // Create dictionary with 5 elements
  inline Dict make_dict(const std::string& n0, const GenericType& x0,
                        const std::string& n1, const GenericType& x1,
                        const std::string& n2, const GenericType& x2,
                        const std::string& n3, const GenericType& x3,
                        const std::string& n4, const GenericType& x4) {
    Dict ret=make_dict(n0, x0, n1, x1, n2, x2, n3, x3);
    ret[n4]=x4;
    return ret;
  }

  // Create dictionary with 6 elements
  inline Dict make_dict(const std::string& n0, const GenericType& x0,
                        const std::string& n1, const GenericType& x1,
                        const std::string& n2, const GenericType& x2,
                        const std::string& n3, const GenericType& x3,
                        const std::string& n4, const GenericType& x4,
                        const std::string& n5, const GenericType& x5) {
    Dict ret=make_dict(n0, x0, n1, x1, n2, x2, n3, x3, n4, x4);
    ret[n5]=x5;
    return ret;
  }

#endif // SWIG


} // namespace casadi


#endif // CASADI_GENERIC_TYPE_HPP
