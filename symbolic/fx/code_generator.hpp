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

#ifndef CODE_GENERATOR_HPP
#define CODE_GENERATOR_HPP

#include "fx.hpp"
#include <sstream>
#include <map>
#include <set>

namespace CasADi{
  
  class CodeGenerator{
  public:

    /// Add an include file optionally using a relative path "..." instead of an absolute path <...>
    void addInclude(const std::string& new_include, bool relative_path = false);

    /// Add an include file optionally using a relative path "..." instead of an absolute path <...>
    int addSparsity(const CRSSparsity& sp);

    /** \brief Get the index of an existing sparsity pattern */
    int getSparsity(const CRSSparsity& sp) const;

    /** \brief Add a dependent function */
    int addDependent(const FX& f);

    /// Flush generated file to a stream
    void flush(std::ostream& s);
    
    /** \brief COPY sparse: y <- x, (see CRSSparsity::set) */
    void generateCopySparse();

    /** Convert in integer to a string */
    static std::string numToString(int n);

    /** \brief  Print to a c file */
    static void printVector(std::ostream &cfile, const std::string& name, const std::vector<int>& v);
  
    /** \brief Print a sparsity pattern to stream */
    static int printSparsity(std::ostream &stream, const CRSSparsity& sp, std::map<const void*,int>& sparsity_index);

    /** \brief Find an existing sparsity pattern */
    static int findSparsity(const CRSSparsity& sp, const std::map<const void*,int>& sparsity_index);

    /** \brief Find an existing dependent function */
    static int findDependent(const FX& f, const std::map<const void*,int>& dependent_index);

    //  private:
    
    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes_;
    std::stringstream auxiliaries_;
    std::stringstream sparsities_;
    std::stringstream dependents_;
    std::stringstream function_;
    std::stringstream finalization_;
    
    // Set of already included header files
    typedef std::set<std::string> StringSet;
    typedef std::map<const void*,int> PointerMap;
    StringSet added_includes_;
    StringSet added_auxiliaries_;
    PointerMap added_sparsities_;
    PointerMap added_dependents_;

  };
  
  
} // namespace CasADi

#endif // CODE_GENERATOR_HPP

