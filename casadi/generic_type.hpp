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

#ifndef GENERIC_TYPE_HPP
#define GENERIC_TYPE_HPP

#include "shared_object.hpp"
#include <string>
#include <vector>

namespace CasADi{

  /** \brief  Types of options */
  enum opt_type { OT_BOOLEAN, OT_INTEGER, OT_REAL, OT_STRING, OT_INTEGERVECTOR, OT_REALVECTOR };

  class GenericTypeInternal;
  
  /** \brief Generic data type
  \author Joel Andersson 
  \date 2010
  Return type when getting an option, can be converted into bool, int, string, vector, etc */
  class GenericType : public SharedObject{
    public:
    GenericType();
    GenericType(bool b);
    GenericType(int i);
    GenericType(double d);
    GenericType(const std::vector<bool>& iv);
    GenericType(const std::vector<int>& iv);
    GenericType(const std::vector<double>& dv);
    GenericType(const std::string& s);
    GenericType(const char s[]);

    /// Implicit typecasting
    #ifndef SWIG
    operator bool() const{ return toBool();} 
    operator int() const{ return toInt();} 
    operator double() const{ return toDouble();}
    operator const std::string& () const{ return toString();}
    operator const std::vector<int>& () const{ return toIntVector();}
    operator const std::vector<double>& () const{ return toDoubleVector();}
    #endif // SWIG
    
    //! \brief Is boolean?
    bool isBool() const;

    //! \brief Is an integer?
    bool isInt() const;
    
    //! \brief Is a double?
    bool isDouble() const;
    
    //! \brief Is a string?
    bool isString() const;

    //! \brief Is a vector of ints?
    bool isIntVector() const;
    
    //! \brief Is a vector of doubles?
    bool isDoubleVector() const;

    //! \brief Convert to boolean
    bool toBool() const;

    //! \brief Convert to int
    int toInt() const;
    
    //! \brief Convert to double
    double toDouble() const;
    
    //! \brief Convert to string
    const std::string& toString() const;

    //! \brief Convert to vector of ints
    const std::vector<int>& toIntVector() const;
    
    //! \brief Convert to vector of doubles
    const std::vector<double>& toDoubleVector() const;

    //! \brief Equality
    bool operator==(const GenericType& op2) const;
    bool operator!=(const GenericType& op2) const;
    
    //! \brief Print
    friend std::ostream& operator<<(std::ostream &stream, const GenericType& ref);
    
    //! \brief Access a member function or object
    //! A regular user is not supposed to use this method.
    GenericTypeInternal* operator->();
    //! \brief Access a member function or object
    //! A regular user is not supposed to use this method.
    const GenericTypeInternal* operator->() const;
    
    /// Check if it is of a certain type
    #ifndef SWIG
    template<typename T>
    bool is_a() const{
      return dynamic_cast<const T*>(get()) != 0;
    }
    #endif // SWIG
    
  };
  

} // namespace CasADi


#endif // GENERIC_TYPE_HPP
