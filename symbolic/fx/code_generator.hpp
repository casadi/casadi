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
    int addDependency(const FX& f);

    /** \brief Get the index of an existing dependency */
    int getDependency(const FX& f) const;

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

    /** \brief  Print to a c file */
    static void printVector(std::ostream &s, const std::string& name, const std::vector<int>& v);

  private:

    /// COPY: y <-x
    void auxCopy();

    /// SWAP: x <-> y
    void auxSwap();
    
    // SCAL: x <- alpha*x
    void auxScal();

    // AXPY: y <- a*x + y
    void auxAxpy();

    // DOT: inner_prod(x,y) -> return
    void auxDot();

    // ASUM: ||x||_1 -> return
    void auxAsum();

    // IAMAX: index corresponding to the entry with the largest absolute value 
    void auxIamax();

    // NRM2: ||x||_2 -> return
    void auxNrm2();

    // FILL: x <- alpha
    void auxFill();
    
    // Sparse matrix-matrix multiplication, the second argument is transposed: z <- z + x*y'
    void auxMmNtSparse();

    /// SIGN
    void auxSign();

    /// COPY sparse: y <- x
    void auxCopySparse();

    /// TRANS: y <- trans(x)
    void auxTrans();

    //  private:
  public:
    
    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes_;
    std::stringstream auxiliaries_;
    std::stringstream sparsities_;
    std::stringstream dependencies_;
    std::stringstream function_;
    std::stringstream finalization_;
    
    // Set of already included header files
    typedef std::map<const void*,int> PointerMap;
    std::set<std::string> added_includes_;
    std::set<Auxiliary> added_auxiliaries_;
    PointerMap added_sparsities_;
    PointerMap added_dependencies_;

  };
  
  
} // namespace CasADi

#endif // CODE_GENERATOR_HPP

