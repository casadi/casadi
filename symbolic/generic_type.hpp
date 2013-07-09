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
#include "casadi_types.hpp"
#include <string>
#include <vector>

namespace CasADi{

  /** \brief  Types of options */
  enum opt_type { OT_BOOLEAN, OT_INTEGER, OT_REAL, OT_STRING, OT_INTEGERVECTOR, OT_BOOLVECTOR, OT_REALVECTOR, OT_STRINGVECTOR, OT_DICTIONARY, OT_NLPSOLVER, OT_LPSOLVER, OT_LINEARSOLVER, OT_INTEGRATOR, OT_QPSOLVER, OT_SOCPSOLVER, OT_QCQPSOLVER, OT_SDPSOLVER,  OT_SDQPSOLVER, OT_IMPLICITFUNCTION, OT_JACOBIANGENERATOR, OT_SPARSITYGENERATOR, OT_FX, OT_VOIDPTR, OT_UNKNOWN};
  
  /** \brief Generic data type
  \author Joel Andersson 
  \date 2010
  Return type when getting an option, can be converted into bool, int, string, vector, etc */
  class GenericType : public SharedObject{
    public:
    GenericType();
    #ifndef SWIG
    GenericType(bool b);
    GenericType(int i);
    GenericType(double d);
    GenericType(const std::string& s);
    GenericType(const std::vector<bool>& iv);
    GenericType(const std::vector<int>& iv);
    GenericType(const std::vector<double>& dv);
    GenericType(const std::vector<std::string>& sv);
    GenericType(const char s[]);
    GenericType(const FX& f);
    #endif // SWIG
    GenericType(const SharedObject& obj);
    //GenericType(const GenericType& obj);    
    
    #ifndef SWIG
    GenericType(void* ptr);
    #endif // SWIG

    /// Get a description of a type
    static std::string get_type_description(const opt_type &type);
    
    /// Get a description of the object's type
    std::string get_description() const { return get_type_description(type_); }

    /// Construct a GenericType given an opt_type
    static GenericType from_type(opt_type type);

    typedef std::map<std::string, GenericType> Dictionary;
    GenericType(const Dictionary& dict);

    /// Creator functions
    GenericType(NLPSolverCreator ptr);
    GenericType(linearSolverCreator ptr);
    GenericType(LPSolverCreator ptr);
    GenericType(integratorCreator ptr);
    GenericType(QPSolverCreator ptr);
    GenericType(SOCPSolverCreator ptr);
    GenericType(SDPSolverCreator ptr);
    GenericType(SDQPSolverCreator ptr);
    GenericType(QCQPSolverCreator ptr);
    GenericType(implicitFunctionCreator ptr);
    GenericType(JacobianGenerator ptr);
    GenericType(SparsityGenerator ptr);
    
    /// Implicit typecasting
    #ifndef SWIG
    operator bool() const{ return toBool();} 
    operator int() const{ return toInt();} 
    operator double() const{ return toDouble();}
    operator const std::string& () const{ return toString();}
    operator const std::vector<int>& () const{ return toIntVector();}
    operator const std::vector<double>& () const{ return toDoubleVector();}
    operator const std::vector<std::string>& () const{ return toStringVector();}
    operator const FX& () const{ return toFX();}
    //operator void*() const;
    operator const std::map<std::string, GenericType>& () const;
    
    operator NLPSolverCreator() const;
    operator linearSolverCreator() const;
    operator LPSolverCreator() const;
    operator integratorCreator() const;
    operator QPSolverCreator() const;
    operator SDPSolverCreator() const;
    operator SDQPSolverCreator() const;
    operator SOCPSolverCreator() const;
    operator QCQPSolverCreator() const;
    operator implicitFunctionCreator() const;
    operator JacobianGenerator() const;
    operator SparsityGenerator() const;
    #endif // SWIG
    
    opt_type getType() const;
    
    bool can_cast_to(opt_type other) const;
    bool can_cast_to(const GenericType& other) const { return can_cast_to(other.type_) ;}
    
    //! \brief Is boolean?
    bool isBool() const;

    //! \brief Is an integer?
    bool isInt() const;
    
    //! \brief Is a double?
    bool isDouble() const;
    
    //! \brief Is a string?
    bool isString() const;

    //! \brief Is an empty vector?
    bool isEmptyVector() const;
    
    //! \brief Is a vector of ints?
    bool isIntVector() const;
    
    //! \brief Is a vector of doubles?
    bool isDoubleVector() const;

    //! \brief Is a vector of strings
    bool isStringVector() const;

    //! \brief Is a shared object?
    bool isSharedObject() const;

    //! \brief Is a shared object?
    bool isDictionary() const;
    
    //! \brief Is a shared object?
    bool isFX() const;
    
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
    
    //! \brief Convert to vector of strings
    const std::vector<std::string>& toStringVector() const;

    //! \brief Convert to shared object
    const SharedObject& toSharedObject() const;

    //! \brief Convert to Dictionary
    const Dictionary& toDictionary() const;

    //! \brief Convert to shared object
    const FX& toFX() const;
    
    //! \brief Convert to void pointer
    void * toVoidPointer() const;

    //! \brief Equality
    bool operator==(const GenericType& op2) const;
    bool operator!=(const GenericType& op2) const;
   
    #ifndef SWIG 
    //! \brief Print
    friend std::ostream& operator<<(std::ostream &stream, const GenericType& ref);
    #endif
        
    /// Check if it is of a certain type (implementation in generic_type_internal.hpp)
    #ifndef SWIG
    template<typename T>
    bool is_a() const;
    #endif // SWIG
    
   private:  
    opt_type type_;
      
    
  };
  

} // namespace CasADi


#endif // GENERIC_TYPE_HPP
