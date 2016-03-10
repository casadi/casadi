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
    CodeGenerator(const Dict& opts = Dict());

    /// Add a function (name generated)
    void add(const Function& f);

#ifndef SWIG
    /// Generate the code to a stream
    void generate(std::ostream& s) const;
#endif // SWIG

    /// Generate a file
    void generate(const std::string& name) const;

    /// Generate a file, return code as string
    std::string generate() const;

    /// Compile and load function
    std::string compile(const std::string& name, const std::string& compiler="gcc -fPIC -O2");

    /// Add an include file optionally using a relative path "..." instead of an absolute path <...>
    void addInclude(const std::string& new_include, bool relative_path=false,
                    const std::string& use_ifdef=std::string());

#ifndef SWIG
    /// Add an external function declaration
    void addExternal(const std::string& new_external);

    // Add a sparsity pattern
    std::string sparsity(const Sparsity& sp);

    // Add a sparsity pattern, get index
    int addSparsity(const Sparsity& sp);

    /** \brief Get the index of an existing sparsity pattern */
    int get_sparsity(const Sparsity& sp) const;

    /** \brief Get or add a constant */
    int getConstant(const std::vector<double>& v, bool allow_adding=false);

    /** \brief Get or add an integer constant */
    int getConstant(const std::vector<int>& v, bool allow_adding=false);

    /** \brief Use simplified signature */
    static bool simplifiedCall(const Function& f);

    /** \brief Generate a call to a function (generic signature) */
    std::string operator()(const Function& f, const std::string& arg,
                           const std::string& res, const std::string& iw,
                           const std::string& w, const std::string& mem="0") const;

    /** \brief Generate a call to a function (simplified signature) */
    std::string operator()(const Function& f, const std::string& arg, const std::string& res) const;

    /** \brief Print a constant in a lossless but compact manner */
    static std::string constant(double v);

    /** \brief Codegen inner product */
    std::string dot(int n, const std::string& x, const std::string& y);

    /** \brief Codegen sparse matrix-matrix multiplication */
    std::string mtimes(const std::string& x, const Sparsity& sp_x,
                       const std::string& y, const Sparsity& sp_y,
                       const std::string& z, const Sparsity& sp_z,
                       const std::string& w, bool tr);

    /** \brief Codegen bilinear form */
    std::string bilin(const std::string& A, const Sparsity& sp_A,
                      const std::string& x, const std::string& y);

    /** \brief Rank-1 update */
    std::string rank1(const std::string& A, const Sparsity& sp_A, const std::string& alpha,
                      const std::string& x, const std::string& y);

    /** \brief Declare a function */
    std::string declare(std::string s);

    /** \brief Auxiliary functions */
    enum Auxiliary {
      AUX_COPY,
      AUX_SWAP,
      AUX_SCAL,
      AUX_AXPY,
      AUX_DOT,
      AUX_BILIN,
      AUX_RANK1,
      AUX_NRM2,
      AUX_IAMAX,
      AUX_FILL,
      AUX_ASUM,
      AUX_SQ,
      AUX_SIGN,
      AUX_MTIMES,
      AUX_PROJECT,
      AUX_TRANS,
      AUX_TO_MEX,
      AUX_FROM_MEX
    };

    /** \brief Add a built-in auxiliary function */
    void addAuxiliary(Auxiliary f);

    /** Convert in integer to a string */
    static std::string to_string(int n);

    /** Get work vector name from index */
    std::string work(int n, int sz) const;

    /** Get work vector element from index */
    std::string workel(int n) const;

    /** \brief  Print int vector to a c file */
    static void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<int>& v);

    /** \brief  Print real vector to a c file */
    static void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<double>& v);

    /** \brief Create a copy operation */
    std::string copy(const std::string& arg, std::size_t n, const std::string& res);

    /** \brief Create a fill operation */
    std::string fill(const std::string& res, std::size_t n, const std::string& v);

    /** \brief Sparse assignment */
    std::string project(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, const Sparsity& sp_res,
                        const std::string& w);

    /** \brief Create matrix in MATLAB's MEX format */
    std::string to_mex(const Sparsity& sp, const std::string& arg);

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

    // Generate header file?
    bool with_header;

    // Are we creating a MEX file?
    bool mex;

    // Verbose codegen?
    bool verbose;

    // Are we generating C++?
    bool cpp;

    // Should we generate a main (allowing evaluation from command line)
    bool main;

    /** \brief Codegen scalar
     * Use the work vector for storing work vector elements of length 1
     * (typically scalar) instead of using local variables
     */
    bool codegen_scalars;

    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes;
    std::stringstream auxiliaries;
    std::stringstream body;
    std::stringstream header;

    // Names of exposed functions
    std::vector<std::string> exposed_fname;

    // Set of already included header files
    typedef std::map<const void*, int> PointerMap;
    std::set<std::string> added_includes_;
    std::set<std::string> added_externals_;
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
#endif // SWIG
  };


} // namespace casadi

#endif // CASADI_CODE_GENERATOR_HPP

