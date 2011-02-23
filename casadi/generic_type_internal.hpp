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

namespace CasADi{
  
  class GenericTypeInternal : public SharedObjectNode{
    public:
    explicit GenericTypeInternal(const std::vector<int>& i_vec);
    explicit GenericTypeInternal(const std::vector<double>& d_vec);
    explicit GenericTypeInternal(const std::string& s);
        
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
    
    //! \brief Length
    int n;
    
    //! \brief Boolean denoting if option type is a string
    bool is_string;
    
    std::vector<int> i_vec;
    std::vector<double> d_vec;
    std::string str;
  };
      

} // namespace CasADi


#endif // GENERIC_TYPE_INTERNAL_HPP
