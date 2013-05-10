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

    /** \brief Get or add a constant */
    int getConstant(const std::vector<double>& v, bool allow_adding=false);

    /** \brief Get or add am integer constant */
    int getConstant(const std::vector<int>& v, bool allow_adding=false);

    /** \brief Add a dependent function */
    int addDependency(const FX& f);

    /** \brief Get the index of an existing dependency */
    int getDependency(const FX& f) const;
    
    /** \brief Print a constant in a lossless but compact manner */
    static void printConstant(std::ostream& s, double v);

    /** \bried Codegen casadi_dot */
    std::string casadi_dot(int n, const std::string& x, int inc_x, const std::string& y, int inc_y);

    /** \brief Auxiliary functions */
    enum Auxiliary{
      // BLAS Level 1
      AUX_COPY,
      AUX_SWAP,
      AUX_SCAL,
      AUX_AXPY,
      AUX_DOT,
      AUX_NRM2,
      AUX_IAMAX,
      AUX_FILL,
      AUX_ASUM,

      // Misc
      AUX_SQ,
      AUX_SIGN,
      AUX_MM_NT_SPARSE,
      AUX_COPY_SPARSE,
      AUX_TRANS
    };
    
    /** \brief Add a built-in axiliary function */
    void addAuxiliary(Auxiliary f);

    /// Flush generated file to a stream
    void flush(std::ostream& s) const;
    
    /** Convert in integer to a string */
    static std::string numToString(int n);

    /** \brief  Print int vector to a c file */
    static void printVector(std::ostream &s, const std::string& name, const std::vector<int>& v);

    /** \brief  Print real vector to a c file */
    static void printVector(std::ostream &s, const std::string& name, const std::vector<double>& v);

    /** \brief Copy a vector to another */
    void copyVector(std::ostream &s, const std::string& arg, std::size_t n, const std::string& res, const std::string& it="i", bool only_if_exists=false) const;

  private:

    /// SQUARE
    void auxSq();

    /// SIGN
    void auxSign();

    //  private:
  public:
    
    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes_;
    std::stringstream auxiliaries_;
    std::stringstream dependencies_;
    std::stringstream function_;
    std::stringstream finalization_;
    
    // Set of already included header files
    typedef std::map<const void*,int> PointerMap;
    std::set<std::string> added_includes_;
    std::set<Auxiliary> added_auxiliaries_;
    PointerMap added_sparsities_;
    PointerMap added_dependencies_;
    std::multimap<size_t,size_t> added_double_constants_;
    std::multimap<size_t,size_t> added_integer_constants_;

    // Constants
    std::vector<std::vector<double> > double_constants_;
    std::vector<std::vector<int> > integer_constants_;

    // Hash a vector
    static size_t hash(const std::vector<double>& v);
    static size_t hash(const std::vector<int>& v);

    // Compare two vectors
    template<typename T>
    static bool equal(const std::vector<T>& v1, const std::vector<T>& v2){
      if(v1.size()!=v2.size()) return false;
      for(int j=0; j<v1.size(); ++j){
        if(v1[j]!=v2[j]) return false;
      }
      return true;
    }
  };
  
  
} // namespace CasADi

#endif // CODE_GENERATOR_HPP

