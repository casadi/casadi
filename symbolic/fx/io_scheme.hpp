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

#ifndef IO_SCHEME_HPP
#define IO_SCHEME_HPP

#include "../shared_object.hpp"
#include "schemes_metadata.hpp"
#include <string>
#include <vector>
#include <iostream>

namespace CasADi{

  // Forward declaration
  class IOSchemeInternal;
  
  /** \brief Class with mapping between names and indices
   *
   * \author Joris Gillis 
   * \date 2013
   */
class IOScheme : public SharedObject{
  
  public:
  
    /// Default constructor
    IOScheme();
    
    /// Constructor with enum
    IOScheme(InputOutputScheme scheme);

    /// Constructor with entry names
    IOScheme(const std::vector<std::string> &entries);
    
    /** \brief  Access functions of the node */
    IOSchemeInternal* operator->();
    const IOSchemeInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
    
    /// Name of the scheme
    std::string name() const;
      
    /// List available entries
    std::string entryNames() const;
      
    /// Get index by entry name
    int index(const std::string &name) const;
      
    /// Number of entries
    int size() const;
    
    /// Get the entry name by index
    std::string entry(int i) const;

    /** \brief Get the entry label by index
    * If scheme is unknown, returns the index as a string
    */
    std::string entryLabel(int i) const;
    
    /// Get the entry enum name by index
    std::string entryEnum(int i) const;
    
    /// Describe the index as an input
    std::string describeInput(int i) const;

    /// Describe the index as an output
    std::string describeOutput(int i) const;
    
    /// Describe the index
    std::string describe(int i) const;
    
    /// Check wether the scheme is known
    bool known() const;
    
    /// Check wether this scheme is compatible with the given size
    int compatibleSize(int size) const;
    
    #ifndef SWIG
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;

    #endif // SWIG  
    
};
  
} // namespace CasADi

#endif // IO_SCHEME_HPP
