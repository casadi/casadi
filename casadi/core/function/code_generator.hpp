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


#ifndef CASADI_CODE_GENERATOR_HPP
#define CASADI_CODE_GENERATOR_HPP

#include "function.hpp"
#include <sstream>
#include <map>
#include <set>

namespace casadi {

  class CASADI_EXPORT CodeGenerator {
  public:
    /// Constructor
    CodeGenerator(const Dictionary& opts = Dictionary());

    /// Add a function
    void addFunction(const Function& f, const std::string& fname);

    /// Generate the code to a stream
    void generate(std::ostream& s) const;

    /// Generate a file
    void generate(const std::string& fname) const;

    /// Add an include file optionally using a relative path "..." instead of an absolute path <...>
    void addInclude(const std::string& new_include, bool relative_path = false);

    // Add a sparsity pattern
    std::string sparsity(const Sparsity& sp);

    // Add a sparsity pattern, get index
    int addSparsity(const Sparsity& sp);

    /** \brief Get the index of an existing sparsity pattern */
    int getSparsity(const Sparsity& sp) const;

    /** \brief Get or add a constant */
    int getConstant(const std::vector<double>& v, bool allow_adding=false);

    /** \brief Get or add an integer constant */
    int getConstant(const std::vector<int>& v, bool allow_adding=false);

    /** \brief Add a dependent function */
    int addDependency(const Function& f);

    /** \brief Get the index of an existing dependency */
    int getDependency(const Function& f) const;

    /** \brief Print a constant in a lossless but compact manner */
    static std::string constant(double v);

    /** \brief Codegen casadi_dot */
    std::string casadi_dot(int n, const std::string& x, int inc_x, const std::string& y, int inc_y);

    /** \brief Auxiliary functions */
    enum Auxiliary {
      AUX_COPY_N,
      AUX_SWAP,
      AUX_SCAL,
      AUX_AXPY,
      AUX_DOT,
      AUX_NRM2,
      AUX_IAMAX,
      AUX_FILL_N,
      AUX_ASUM,
      AUX_SQ,
      AUX_SIGN,
      AUX_MM_SPARSE,
      AUX_PROJECT,
      AUX_TRANS,
      AUX_TO_MEX,
      AUX_FROM_MEX
    };

    /** \brief Add a built-in auxiliary function */
    void addAuxiliary(Auxiliary f);

    /** Convert in integer to a string */
    static std::string numToString(int n);

    /** Get work vector name from index */
    static std::string work(int n);

    /** Get work vector element from index */
    static std::string workelement(int n);

    /** \brief  Print int vector to a c file */
    static void printVector(std::ostream &s, const std::string& name, const std::vector<int>& v);

    /** \brief  Print real vector to a c file */
    static void printVector(std::ostream &s, const std::string& name, const std::vector<double>& v);

    /** \brief Create a copy_n operation */
    std::string copy_n(const std::string& arg, std::size_t arg_off, std::size_t n,
                       const std::string& res, std::size_t res_off);

    /** \brief Create a fill_n operation */
    std::string fill_n(const std::string& res, std::size_t res_off, std::size_t n,
                       const std::string& v);

    /** \brief Sparse assignment */
    std::string project(const std::string& arg, std::size_t arg_off, const Sparsity& sp_arg,
                        const std::string& res, std::size_t res_off, const Sparsity& sp_res,
                        const std::string& w);

    /** \brief Create matrix in MATLAB's MEX format */
    std::string to_mex(const Sparsity& sp, const std::string& data="0");

    /** \brief Get matrix from MATLAB's MEX format */
    std::string from_mex(std::string& arg,
                         const std::string& res, std::size_t res_off, const Sparsity& sp_res,
                         const std::string& w);

    /** \brief Assignment */
    static void assign(std::ostream &s, const std::string& lhs, const std::string& rhs);

    /** \brief Printf */
    std::string printf(const std::string& str,
                       const std::vector<std::string>& arg=std::vector<std::string>());
    std::string printf(const std::string& str, const std::string& arg1);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2,
                       const std::string& arg3);
  private:

    /// SQUARE
    void auxSq();

    /// SIGN
    void auxSign();

    //  private:
  public:

    /// \cond INTERNAL

    // Real-type used for the codegen
    std::string real_t;

    // Are we creating a MEX file?
    bool mex;

    // C++ guards (making the code valid C++)
    bool cpp_guards;

    // Should we generate a main (allowing evaluation from command line)
    bool main;

    // Prefix
    std::string prefix;

    // Include file
    std::string include;

    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes;
    std::stringstream auxiliaries;
    std::stringstream functions;

    // Set of already included header files
    typedef std::map<const void*, int> PointerMap;
    std::set<std::string> added_includes_;
    std::set<Auxiliary> added_auxiliaries_;
    PointerMap added_sparsities_;
    PointerMap added_dependencies_;
    std::multimap<size_t, size_t> added_double_constants_;
    std::multimap<size_t, size_t> added_integer_constants_;

    // Constants
    std::vector<std::vector<double> > double_constants_;
    std::vector<std::vector<int> > integer_constants_;

    // Hash a vector
    static size_t hash(const std::vector<double>& v);
    static size_t hash(const std::vector<int>& v);

    // Compare two vectors
    template<typename T>
    static bool equal(const std::vector<T>& v1, const std::vector<T>& v2) {
      if (v1.size()!=v2.size()) return false;
      for (int j=0; j<v1.size(); ++j) {
        if (v1[j]!=v2[j]) return false;
      }
      return true;
    }
    /// \endcond
  };


} // namespace casadi

#endif // CASADI_CODE_GENERATOR_HPP

