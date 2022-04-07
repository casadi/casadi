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

#include <map>
#include <set>
#include <sstream>

namespace casadi {

  /** \brief Helper class for C code generation
      \author Joel Andersson
      \date 2016
  */
  class CASADI_EXPORT CodeGenerator {
  public:
    /// Constructor
    CodeGenerator(const std::string& name, const Dict& opts = Dict());

    /// Add a function (name generated)
    void add(const Function& f, bool with_jac_sparsity=false);

#ifndef SWIG
    /// Generate the code to a stream
    void dump(std::ostream& s);
#endif // SWIG

    /// Generate a file, return code as string
    std::string dump();

    /** \brief Generate file(s)
      The "prefix" argument will be prepended to the generated files and may
      be a directory or a file prefix.
      returns the filename
    */
    std::string generate(const std::string& prefix="");

    /// Add an include file optionally using a relative path "..." instead of an absolute path <...>
    void add_include(const std::string& new_include, bool relative_path=false,
                    const std::string& use_ifdef=std::string());

#ifndef SWIG
    /// Add a function dependency
    std::string add_dependency(const Function& f);

    /// Add an external function declaration
    void add_external(const std::string& new_external);

    /// Get a shorthand
    std::string shorthand(const std::string& name) const;

    /// Add/get a shorthand
    std::string shorthand(const std::string& name, bool allow_adding=true);

    // Add a sparsity pattern
    std::string sparsity(const Sparsity& sp);

    // Add a sparsity pattern, get index
    casadi_int add_sparsity(const Sparsity& sp);

    /** \brief Get the index of an existing sparsity pattern */
    casadi_int get_sparsity(const Sparsity& sp) const;

    /** \brief Get or add a constant */
    casadi_int get_constant(const std::vector<double>& v, bool allow_adding=false);

    /** \brief Get or add an integer constant */
    casadi_int get_constant(const std::vector<casadi_int>& v, bool allow_adding=false);

    /** \brief Represent an array constant; adding it when new */
    std::string constant(const std::vector<casadi_int>& v);

    /** \brief Represent an array constant; adding it when new */
    void constant_copy(const std::string& var_name, const std::vector<casadi_int>& v, const std::string& type="casadi_int");

    /** \brief Represent an array constant; adding it when new */
    std::string constant(const std::vector<double>& v);

    /** \brief Allocate file scope double read-only memory */
    void define_rom_double(const void* id, casadi_int size);

    /** \brief Access file scope double read-only memory */
    std::string rom_double(const void* id) const;

    /** \brief Allocate file scope integer read-only memory */
    void define_rom_integer(const void* id, casadi_int size);

    /** \brief Access file scope integer read-only memory */
    std::string rom_integer(const void* id) const;

    /** \brief Generate a call to a function (generic signature) */
    std::string operator()(const Function& f, const std::string& arg,
                           const std::string& res, const std::string& iw,
                           const std::string& w);

    /** \brief Print a string to buffer  */
    CodeGenerator& operator<<(const std::string& s);

    /** \brief Print without newline characters */
    void print_formatted(const std::string& s);

    /** \brief Print an arbitrary type to buffer */
    template<typename T>
    CodeGenerator& operator<<(T s) {
      std::stringstream ss;
      ss << s;
      return (*this) << ss.str();
    }

    /** \brief Flush the buffer to a stream of choice */
    void flush(std::ostream &s);

    /** \brief Declare a local variable */
    void local(const std::string& name, const std::string& type, const std::string& ref="");

    /** \brief Enter a local scope */
    void scope_enter();

    /** \brief Exit a local scope */
    void scope_exit();

    /** \brief Declare a work vector element */
    std::string sx_work(casadi_int i);

    /** \brief Specify the default value for a local variable */
    void init_local(const std::string& name, const std::string& def);

    /** \brief Increase indentation */
    void indent() {current_indent_++;}

    /** \brief Decrease indentation */
    void unindent() {current_indent_--;}

    /** \brief Avoid stack? */
    bool avoid_stack() { return avoid_stack_;}

    /** \brief Print a constant in a lossless but compact manner */
    std::string constant(double v);
    std::string constant(casadi_int v);

    /** \brief Print an intializer */
    std::string initializer(const std::vector<double>& v);
    std::string initializer(const std::vector<casadi_int>& v);

    /** \brief Sanitize source files for codegen */
    std::string sanitize_source(const std::string& src,
                                const std::vector<std::string>& inst,
                                bool add_shorthand=true);

    /** \brief Codegen inner product */
    std::string dot(casadi_int n, const std::string& x, const std::string& y);

    /** \brief Codegen sparse matrix-vector multiplication */
    std::string mv(const std::string& x, const Sparsity& sp_x,
                   const std::string& y, const std::string& z, bool tr);

    /** \brief Codegen dense matrix-vector multiplication */
    std::string mv(const std::string& x, casadi_int nrow_x, casadi_int ncol_x,
                   const std::string& y, const std::string& z, bool tr);

    /** \brief Codegen axpy: y += a*x */
    std::string axpy(casadi_int n, const std::string& a,
                     const std::string& x, const std::string& y);

    /** \brief Codegen axpy: x *= alpha */
    std::string scal(casadi_int n, const std::string& alpha, const std::string& x);

    /** \brief Codegen sparse matrix-matrix multiplication */
    std::string mtimes(const std::string& x, const Sparsity& sp_x,
                       const std::string& y, const Sparsity& sp_y,
                       const std::string& z, const Sparsity& sp_z,
                       const std::string& w, bool tr);

    /** \brief Codegen lower triangular solve */
    std::string trilsolve(const Sparsity& sp_x, const std::string& x, const std::string& y,
                          bool tr, bool unity, casadi_int nrhs);

    /** \brief Codegen upper triangular solve */
    std::string triusolve(const Sparsity& sp_x, const std::string& x, const std::string& y,
                          bool tr, bool unity, casadi_int nrhs);

    /** \brief Codegen bilinear form */
    std::string bilin(const std::string& A, const Sparsity& sp_A,
                      const std::string& x, const std::string& y);

    /** \brief Rank-1 update */
    std::string rank1(const std::string& A, const Sparsity& sp_A, const std::string& alpha,
                      const std::string& x, const std::string& y);

    /** \brie LogSumExp */
    std::string logsumexp(const std::string& A, casadi_int n);

    /** \brief Multilinear interpolation */
    std::string interpn(const std::string& res, casadi_int ndim, const std::string& grid,
                        const std::string& offset,
                        const std::string& values, const std::string& x,
                        const std::string& lookup_mode, casadi_int m,
                        const std::string& iw, const std::string& w);

    /** \brief Multilinear interpolation - calculate gradient */
    std::string interpn_grad(const std::string& grad,
      casadi_int ndim, const std::string& grid,
      const std::string& offset,
      const std::string& values, const std::string& x,
      const std::string& lookup_mode, casadi_int m,
      const std::string& iw, const std::string& w);

    /** \brief Transpose */
    std::string trans(const std::string& x, const Sparsity& sp_x,
      const std::string& y, const Sparsity& sp_y, const std::string& iw);

    /** \brief QR factorization */
    std::string qr(const std::string& sp, const std::string& A,
                   const std::string& w, const std::string& sp_v,
                   const std::string& v, const std::string& sp_r,
                   const std::string& r, const std::string& beta,
                   const std::string& prinv, const std::string& pc);

    /** \brief QR solve */
    std::string qr_solve(const std::string& x, casadi_int nrhs, bool tr,
                         const std::string& sp_v, const std::string& v,
                         const std::string& sp_r, const std::string& r,
                         const std::string& beta, const std::string& prinv,
                         const std::string& pc, const std::string& w);

    /** \\brief LSQR solve */
    std::string lsqr_solve(const std::string& A, const std::string&x,
                          casadi_int nrhs, bool tr, const std::string& sp, const std::string& w);

    /** \brief LDL factorization */
    std::string ldl(const std::string& sp_a, const std::string& a,
                   const std::string& sp_lt, const std::string& lt,
                   const std::string& d, const std::string& p,
                   const std::string& w);

    /** \brief LDL solve */
    std::string ldl_solve(const std::string& x, casadi_int nrhs,
                         const std::string& sp_lt, const std::string& lt,
                         const std::string& d, const std::string& p,
                         const std::string& w);

    /** \brief fmax */
    std::string fmax(const std::string& x, const std::string& y);

    /** \brief fmin */
    std::string fmin(const std::string& x, const std::string& y);

    /** \brief mmax */
    std::string mmax(const std::string& x, casadi_int n, bool is_dense);

    /** \brief mmin */
    std::string mmin(const std::string& x, casadi_int n, bool is_dense);

    /** \brief vfmax */
    std::string vfmax(const std::string& x, casadi_int n, const std::string& y);

    /** \brief vfmin */
    std::string vfmin(const std::string& x, casadi_int n, const std::string& y);

    /** \brief vfmax */
    std::string vfmax(const std::string& x, const std::string& n, const std::string& y);

    /** \brief vfmin */
    std::string vfmin(const std::string& x, const std::string& n, const std::string& y);

    /** \brief max */
    std::string max(const std::string& x, const std::string& y);

    /** \brief min */
    std::string min(const std::string& x, const std::string& y);

    /** \brief norm_inf */
    std::string norm_inf(casadi_int n, const std::string& x);

    /** \brief max_viol */
    std::string max_viol(casadi_int n, const std::string& x,
      const std::string& lb, const std::string& ub);

    /** \brief sum_viol */
    std::string sum_viol(casadi_int n, const std::string& x,
      const std::string& lb, const std::string& ub);

    /** \brief bound_consistency */
    std::string bound_consistency(casadi_int n, const std::string& x,
      const std::string& lam, const std::string& lbx, const std::string& ubx);

    /** \brief lb_eig */
    std::string lb_eig(const Sparsity& sp_h, const std::string& h);

    /** \brief regularize */
    std::string regularize(const Sparsity& sp_h, const std::string& h, const std::string& reg);

    /** \brief convexify */
    std::string convexify_eval(const ConvexifyData& d,
      const std::string& Hin, const std::string& Hout, const std::string& iw, const std::string& w);

    /** \brief low */
    std::string low(const std::string& x, const std::string& grid,
      casadi_int ng, casadi_int lookup_mode);

    /** \brief Declare a function */
    std::string declare(std::string s);

    /** \brief Write a comment line (ignored if not verbose) */
    void comment(const std::string& s);

    /** \brief Auxiliary functions */
    enum Auxiliary {
      AUX_COPY,
      AUX_CVX,
      AUX_CONVEXIFY,
      AUX_SWAP,
      AUX_SCAL,
      AUX_AXPY,
      AUX_DOT,
      AUX_BILIN,
      AUX_RANK1,
      AUX_NORM_1,
      AUX_NORM_2,
      AUX_NORM_INF,
      AUX_IAMAX,
      AUX_CLEAR,
      AUX_FILL,
      AUX_MV,
      AUX_MV_DENSE,
      AUX_MTIMES,
      AUX_TRILSOLVE,
      AUX_TRIUSOLVE,
      AUX_PROJECT,
      AUX_TRI_PROJECT,
      AUX_DENSIFY,
      AUX_SPARSIFY,
      AUX_TRANS,
      AUX_TO_MEX,
      AUX_FROM_MEX,
      AUX_INTERPN,
      AUX_INTERPN_GRAD,
      AUX_FLIP,
      AUX_INTERPN_WEIGHTS,
      AUX_LOW,
      AUX_INTERPN_INTERPOLATE,
      AUX_DE_BOOR,
      AUX_ND_BOOR_EVAL,
      AUX_FINITE_DIFF,
      AUX_QR,
      AUX_QP,
      AUX_NLP,
      AUX_SQPMETHOD,
      AUX_LDL,
      AUX_NEWTON,
      AUX_TO_DOUBLE,
      AUX_TO_INT,
      AUX_CAST,
      AUX_SQ,
      AUX_SIGN,
      AUX_IF_ELSE,
      AUX_PRINTF,
      AUX_FMIN,
      AUX_FMAX,
      AUX_FABS,
      AUX_MIN,
      AUX_MAX,
      AUX_VFMIN,
      AUX_VFMAX,
      AUX_MAX_VIOL,
      AUX_SUM_VIOL,
      AUX_REGULARIZE,
      AUX_INF,
      AUX_NAN,
      AUX_REAL_MIN,
      AUX_ISINF,
      AUX_BOUNDS_CONSISTENCY,
      AUX_LSQR,
      AUX_FILE_SLURP,
      AUX_CACHE,
      AUX_LOG1P,
      AUX_EXPM1,
      AUX_HYPOT,
      AUX_MMIN,
      AUX_MMAX,
      AUX_LOGSUMEXP
    };

    /** \brief Add a built-in auxiliary function */
    void add_auxiliary(Auxiliary f, const std::vector<std::string>& inst = {"casadi_real"});

    /** \brief Add io sparsity patterns of a function */
    void add_io_sparsities(const std::string& name,
                           const std::vector<Sparsity>& sp_in,
                           const std::vector<Sparsity>& sp_out);

    /** Get work vector name from index */
    std::string work(casadi_int n, casadi_int sz) const;

    /** Get work vector element from index */
    std::string workel(casadi_int n) const;

    /** Declare an array */
    static std::string array(const std::string& type, const std::string& name, casadi_int len,
                             const std::string& def=std::string());

    /** \brief  Print casadi_int vector to a c file */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<casadi_int>& v);

    /** \brief  Print real vector to a c file */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<double>& v);

    /** \brief Create a copy operation */
    std::string copy(const std::string& arg, std::size_t n, const std::string& res);
    void copy_check(const std::string& arg, std::size_t n, const std::string& res,
      bool check_lhs=true, bool check_rhs=true);
    void copy_default(const std::string& arg, std::size_t n, const std::string& res,
      const std::string& def,  bool check_rhs=true);

    /** \brief Create a fill operation */
    std::string fill(const std::string& res, std::size_t n, const std::string& v);

    /** \brief Create a fill operation */
    std::string clear(const std::string& res, std::size_t n);

    /** \brief Refer to argument */
    std::string arg(casadi_int i) const;

    /** \brief Refer to resuly */
    std::string res(casadi_int i) const;

    /** \brief Access thread-local memory */
    std::string mem(const Function& f);

    /** \brief Sparse assignment */
    std::string project(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, const Sparsity& sp_res,
                        const std::string& w);

    /** \brief Project triangular part */
    std::string tri_project(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, bool lower);

    /** \brief Densify */
    std::string densify(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, bool tr=false);

    /** \brief Sparsify */
    std::string sparsify(const std::string& arg, const std::string& res,
                         const Sparsity& sp_res, bool tr=false);

    /** \brief Create matrix in MATLAB's MEX format */
    std::string to_mex(const Sparsity& sp, const std::string& arg);

    /** \brief Get matrix from MATLAB's MEX format */
    std::string from_mex(std::string& arg,
                         const std::string& res, std::size_t res_off, const Sparsity& sp_res,
                         const std::string& w);

    /** \brief Printf */
    std::string printf(const std::string& str,
                       const std::vector<std::string>& arg=std::vector<std::string>());
    std::string printf(const std::string& str, const std::string& arg1);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2,
                       const std::string& arg3);

    /** \brief Print an operation to a c file */
    std::string print_op(Operation op, const std::string& a0);
    std::string print_op(Operation op, const std::string& a0, const std::string& a1);

    /** \brief Slurp a file */
    std::string file_slurp(const std::string& fname, casadi_int n, const std::string& a);

    /** \brief cache check */
    std::string cache_check(const std::string& key, const std::string& cache,
        const std::string& loc, casadi_int stride, casadi_int sz, casadi_int key_sz,
        const std::string& val);
  private:

    /// Print file header
    void file_open(std::ofstream& f, const std::string& name) const;

    /// Print file header
    void file_close(std::ofstream& f) const;

    // Generate casadi_real definition
    void generate_casadi_real(std::ostream &s) const;

    // Generate casadi_int definition
    void generate_casadi_int(std::ostream &s) const;

    // Generate mex entry point
    void generate_mex(std::ostream &s) const;

    // Generate main entry point
    void generate_main(std::ostream &s) const;

    // Generate export symbol macros
    void generate_export_symbol(std::ostream &s) const;

    // Generate import symbol macros
    void generate_import_symbol(std::ostream &s) const;

    //  private:
  public:
    /// \cond INTERNAL

    // Name of generated file
    std::string name, suffix;

    // Real-type used for the codegen
    std::string casadi_real_type;

    // Int-type used for the codegen
    std::string casadi_int_type;

    // Should we create a memory entry point?
    bool with_mem;

    // Generate header file?
    bool with_header;

    // Are we creating a MEX file?
    bool mex;

    // Verbose codegen?
    bool verbose;

    // Verbose runtime?
    bool verbose_runtime;

    // Are we generating C++?
    bool cpp;

    // Should we generate a main (allowing evaluation from command line)
    bool main;

    // Should we include mayth library?
    bool include_math;

    // Do we want to be lean on stack usage?
    bool avoid_stack_;

    std::string infinity, nan, real_min;

    /** \brief Codegen scalar
     * Use the work vector for storing work vector elements of length 1
     * (typically scalar) instead of using local variables
     */
    bool codegen_scalars;

    // Have a flag for exporting/importing symbols
    bool with_export, with_import;

    // Prefix symbols in DLLs?
    std::string dll_export, dll_import;

    // Prefix
    std::string prefix;

    // Stringstreams holding the different parts of the file being generated
    std::stringstream includes;
    std::stringstream auxiliaries;
    std::stringstream body;
    std::stringstream header;
    std::stringstream buffer;

    // Are we at a new line?
    bool newline_;

    // Indentation
    casadi_int indent_;
    casadi_int current_indent_;

    // Names of exposed functions
    std::vector<std::string> exposed_fname;

    // Code generated sparsities
    std::set<std::string> sparsity_meta;

    // Set of already included header files
    std::set<std::string> added_includes_;
    std::set<std::string> added_externals_;
    std::set<std::string> added_shorthands_;
    std::multimap<Auxiliary, std::vector<std::string>> added_auxiliaries_;
    std::multimap<size_t, size_t> added_double_constants_;
    std::multimap<size_t, size_t> added_integer_constants_;
    std::map<std::string, std::pair<std::string, std::string> > local_variables_;
    std::map<std::string, std::string> local_default_;
    std::map<const void *, casadi_int> file_scope_double_;
    std::map<const void *, casadi_int> file_scope_integer_;

    // Added functions
    struct FunctionMeta {
      // The function object
      Function f;
      // Name in codegen
      std::string codegen_name;
    };
    std::vector<FunctionMeta> added_functions_;

    // Constants
    std::vector<std::vector<double> > double_constants_;
    std::vector<std::vector<casadi_int> > integer_constants_;

    // Does any function need thread-local memory?
    bool needs_mem_;

    // Hash a vector
    static size_t hash(const std::vector<double>& v);
    static size_t hash(const std::vector<casadi_int>& v);

    // Compare two vectors
    template<typename T>
    static bool equal(const std::vector<T>& v1, const std::vector<T>& v2) {
      if (v1.size()!=v2.size()) return false;
      for (casadi_int j=0; j<v1.size(); ++j) {
        if (v1[j]!=v2[j]) return false;
      }
      return true;
    }
    /// \endcond
#endif // SWIG
  };


} // namespace casadi

#endif // CASADI_CODE_GENERATOR_HPP
