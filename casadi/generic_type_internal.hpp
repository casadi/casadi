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

#ifndef GENERIC_TYPE_INTERNAL_HPP
#define GENERIC_TYPE_INTERNAL_HPP

#include "generic_type.hpp"
#include "stl_vector_tools.hpp"

namespace CasADi{
  
  class GenericTypeInternal : public SharedObjectNode{
    public:
    explicit GenericTypeInternal(){}
    virtual ~GenericTypeInternal(){}
        
    //! \brief Convert to boolean
    bool toBool() const;
    //! \brief Convert to int
    int toInt() const;
    //! \brief Convert to double
    double toDouble() const;
    //! \brief Convert to string
    virtual const std::string& toString() const{ throw CasadiException("not a string"); }
    //! \brief Convert to vector of ints
    virtual const std::vector<int>& toIntVector() const{ throw CasadiException("not an int vector"); };
    //! \brief Convert to vector of doubles
    virtual const std::vector<double>& toDoubleVector() const{ throw CasadiException("not a double vector"); };

    // Printing
    virtual void print(std::ostream &stream) const = 0;
  };
      
  class StringType : public GenericTypeInternal{
    public:
      explicit StringType(const std::string& d) : d_(d){}
      virtual ~StringType(){}
      virtual const std::string& toString() const{ return d_; }
      virtual void print(std::ostream &stream) const{ stream << d_; }
      std::string d_;
  };
   
  class DoubleVectorType : public GenericTypeInternal{
    public:
      explicit DoubleVectorType(const std::vector<double>& d) : d_(d){}
      virtual ~DoubleVectorType(){}
      virtual const std::vector<double>& toDoubleVector() const{ return d_; }
      virtual void print(std::ostream &stream) const{ stream << d_; }
      std::vector<double> d_;
  };
   
  class IntVectorType : public GenericTypeInternal{
    public:
      explicit IntVectorType(const std::vector<int>& d) : d_(d){}
      virtual ~IntVectorType(){}
      virtual const std::vector<int>& toIntVector() const{ return d_; }
      virtual void print(std::ostream &stream) const{ stream << d_; }
      std::vector<int> d_;
  };
  

} // namespace CasADi


#endif // GENERIC_TYPE_INTERNAL_HPP
