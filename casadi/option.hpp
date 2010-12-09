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

#ifndef OPTION
#define OPTION

#include "shared_object.hpp"
#include <string>
#include <vector>

namespace CasADi{

  /** \brief  Types of options */
  enum opt_type { OT_BOOLEAN, OT_INTEGER, OT_REAL, OT_STRING, OT_INTEGERVECTOR, OT_REALVECTOR };

  class OptionNode;
  
  /** \brief Option data type for the OptionsFunctionality class
  \author Joel Andersson 
  \date 2010
  Return type when getting an option, can be converted into bool, int, string, vector, etc */
  class Option : public SharedObject{
    public:
    Option();
    Option(int i);
    Option(double d);
    Option(const std::vector<bool>& iv);
    Option(const std::vector<int>& iv);
    Option(const std::vector<double>& dv);
    Option(const std::string& s);
    Option(const char s[]);
    
    //! \brief Convert to boolean
    bool toBool() const;
    //! \brief Convert to int
    int toInt() const;
    //! \brief Convert to double
    double toDouble() const;
    //! \brief Convert to string
    std::string toString() const;
    //! \brief Convert to vector of ints
    const std::vector<int>& toIntVector() const;
    //! \brief Convert to vector of doubles
    const std::vector<double>& toDoubleVector() const;

    //! \brief Equality
    friend bool operator==(const Option& op1, const Option& op2);
    friend bool operator!=(const Option& op1, const Option& op2);
    
    //! \brief Print
    friend std::ostream& operator<<(std::ostream &stream, const Option& ref);
    
    //! \brief Access a member function or object
    //! A regular user is not supposed to use this method.
    OptionNode* operator->();
    //! \brief Access a member function or object
    //! A regular user is not supposed to use this method.
    const OptionNode* operator->() const;
  };
  
  /** \brief N class for storing Option data for the OptionsFunctionality class
  \author Joel Andersson 
  \date 2010
  A regular user is not supposed to work with this Node class. Use Option directly.
   */
  class OptionNode : public SharedObjectNode{
    public:
    explicit OptionNode(const std::vector<int>& i_vec);
    explicit OptionNode(const std::vector<double>& d_vec);
    explicit OptionNode(const std::string& s);
        
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


#endif // OPTION
