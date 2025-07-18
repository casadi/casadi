/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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

      \identifier{ru} */
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

        \identifier{rv} */
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

    /* Add a sparsity pattern
    *
    * \param canonical If true, request canonical form,
    * as opposed to potential dense abbreviation
    */
    std::string sparsity(const Sparsity& sp, bool canonical=true);

    /* Add a sparsity pattern, get index
    *
    * \param canonical If true, request canonical form,
    * as opposed to potential dense abbreviation
    */
    casadi_int add_sparsity(const Sparsity& sp, bool canonical=true);

    /** \brief Get the index of an existing sparsity pattern

        \identifier{rw} */
    casadi_int get_sparsity(const Sparsity& sp) const;

    /** \brief Get or add a constant

        \identifier{rx} */
    casadi_int get_constant(const std::vector<double>& v, bool allow_adding=false);

    /** \brief Get or add an integer constant

        \identifier{ry} */
    casadi_int get_constant(const std::vector<casadi_int>& v, bool allow_adding=false);

    /** \brief Get or add a char constant

        \identifier{27l} */
    casadi_int get_constant(const std::vector<char>& v, bool allow_adding=false);

    /** \brief Get or add a vector<string> constant

        \identifier{27z} */
    casadi_int get_constant(const std::vector<std::string>& v, bool allow_adding=false);

    /** \brief Represent an array constant; adding it when new

        \identifier{rz} */
    std::string constant(const std::vector<casadi_int>& v);

    /** \brief Represent an array constant; adding it when new

        \identifier{255} */
    std::string constant(const std::vector<int>& v) {
        return constant(vector_static_cast<casadi_int>(v));
    }

    /** \brief Represent an array constant; adding it when new

        \identifier{s0} */
    void constant_copy(
        const std::string& var_name,
        const std::vector<casadi_int>& v,
        const std::string& type="casadi_int");

    /** \brief Represent an array constant; adding it when new

        \identifier{s1} */
    std::string constant(const std::vector<double>& v);

    /** \brief Represent an array constant; adding it when new

        \identifier{27m} */
    std::string constant(const std::vector<char>& v);

    /** \brief Represent an array constant; adding it when new

        \identifier{280} */
    std::string constant(const std::vector<std::string>& v);

    /** \brief Allocate file scope double read-only memory

        \identifier{s2} */
    void define_rom_double(const void* id, casadi_int size);

    /** \brief Access file scope double read-only memory

        \identifier{s3} */
    std::string rom_double(const void* id) const;

    /** \brief Allocate file scope integer read-only memory

        \identifier{s4} */
    void define_rom_integer(const void* id, casadi_int size);

    /** \brief Access file scope integer read-only memory

        \identifier{s5} */
    std::string rom_integer(const void* id) const;

    /** \brief Allocate file scope double writeable memory

        \identifier{2aw} */
    void define_pool_double(const std::string& name, const std::vector<double>& def);

    /** \brief Access file scope double writeable memory

        \identifier{2ax} */
    std::string pool_double(const std::string& name) const;

    /** \brief Setup a callback

        \identifier{27s} */
    void setup_callback(const std::string& s, const Function& f);

    /** \brief Generate a call to a function (generic signature)

        \identifier{s6} */
    std::string operator()(const Function& f, const std::string& arg,
                           const std::string& res, const std::string& iw,
                           const std::string& w, const std::string& failure_ret="1");

    /** \brief Print a string to buffer

        \identifier{s7} */
    CodeGenerator& operator<<(const std::string& s);

    /** \brief Print without newline characters

        \identifier{s8} */
    void print_formatted(const std::string& s);

    /** \brief Print an arbitrary type to buffer

        \identifier{s9} */
    template<typename T>
    CodeGenerator& operator<<(T s) {
      std::stringstream ss;
      ss << s;
      return (*this) << ss.str();
    }

    /** \brief Flush the buffer to a stream of choice

        \identifier{sa} */
    void flush(std::ostream &s);

    /** \brief Declare a local variable

        \identifier{sb} */
    void local(const std::string& name, const std::string& type, const std::string& ref="");

    /** \brief Enter a local scope

        \identifier{sc} */
    void scope_enter();

    /** \brief Exit a local scope

        \identifier{sd} */
    void scope_exit();

    /** \brief Declare a work vector element

        \identifier{se} */
    std::string sx_work(casadi_int i);

    /** \brief Specify the default value for a local variable

        \identifier{sf} */
    void init_local(const std::string& name, const std::string& def);

    /** \brief Increase indentation

        \identifier{sg} */
    void indent() {current_indent_++;}

    /** \brief Decrease indentation

        \identifier{sh} */
    void unindent() {current_indent_--;}

    /** \brief Avoid stack?

        \identifier{si} */
    bool avoid_stack() const { return avoid_stack_;}

    /** \brief Print a constant in a lossless but compact manner

        \identifier{sj} */
    std::string constant(double v);
    std::string constant(casadi_int v);
    std::string constant(const std::string& v);
    std::string constant(char v);

    std::string format_padded(casadi_int i) const;

    std::string zeros(casadi_int sz);
    std::string ones(casadi_int sz);

    /** \brief Print an initializer

        \identifier{sk} */
    template <typename T>
    std::string initializer(const std::vector<T>& v) {
        std::stringstream s;
        if (v.size() > max_initializer_elements_per_line) {
            s << "\n  ";
        }

        s << "{";
        for (casadi_int i = 0; i < v.size(); ++i) {
            if (i != 0) {
                if (max_initializer_elements_per_line > 1 &&
                    i % max_initializer_elements_per_line == 0) {
                    s << ",\n  ";
                } else {
                    s << ", ";
                }
            }
            s << constant(v[i]);
        }
        s << "}";
        return s.str();
    }

    /** \brief Sanitize source files for codegen

        \identifier{sl} */
    std::string sanitize_source(const std::string& src,
                                const std::vector<std::string>& inst,
                                bool add_shorthand=true);

    /** \brief Codegen inner product

        \identifier{sm} */
    std::string dot(casadi_int n, const std::string& x, const std::string& y);

    /** \brief Codegen sparse matrix-vector multiplication

        \identifier{sn} */
    std::string mv(const std::string& x, const Sparsity& sp_x,
                   const std::string& y, const std::string& z, bool tr);

    /** \brief Codegen dense matrix-vector multiplication

        \identifier{so} */
    std::string mv(const std::string& x, casadi_int nrow_x, casadi_int ncol_x,
                   const std::string& y, const std::string& z, bool tr);

    /** \brief Codegen axpy: y += a*x

        \identifier{sp} */
    std::string axpy(casadi_int n, const std::string& a,
                     const std::string& x, const std::string& y);

    /**
     * @brief Codegen clip_min: Clips the smaller entries in a vector than min
     * to the min
     
     */
    std::string clip_min(const std::string& x, casadi_int n,
                                  const std::string& min, const std::string& mask);

    /**
     * @brief Codegen clip_max: Clips the larger entries in a vector than max
     * to the max
     
     */
    std::string clip_max(const std::string& x, casadi_int n,
                                    const std::string& min, const std::string& mask);

    /**
     * @brief Codegen vector_fmax: Takes vectorwise max of a vector and writes
     * the result to second vector
     */
    std::string vector_fmax(casadi_int n, const std::string& x,
                                    const std::string& y, const std::string& z);

    /**
     * @brief Codegen vector_fmin: Takes vectorwise min of a vector and writes
     * the result to second vector
     */
    std::string vector_fmin(casadi_int n, const std::string& x,
                                    const std::string& y, const std::string& z);

    /**
     * @brief codegen masked_norm_inf: The mask tells what entry is used in the
     * inf-norm.
     */
    std::string masked_norm_inf(casadi_int n, const std::string& x,
                                    const std::string& mask);


    /** \brief What does scal do??

        \identifier{sq} */
    std::string scal(casadi_int n, const std::string& alpha, const std::string& x);

    /** \brief Codegen sparse matrix-matrix multiplication

        \identifier{sr} */
    std::string mtimes(const std::string& x, const Sparsity& sp_x,
                       const std::string& y, const Sparsity& sp_y,
                       const std::string& z, const Sparsity& sp_z,
                       const std::string& w, bool tr);

    /** \brief Codegen lower triangular solve

        \identifier{ss} */
    std::string trilsolve(const Sparsity& sp_x, const std::string& x, const std::string& y,
                          bool tr, bool unity, casadi_int nrhs);

    /** \brief Codegen upper triangular solve

        \identifier{st} */
    std::string triusolve(const Sparsity& sp_x, const std::string& x, const std::string& y,
                          bool tr, bool unity, casadi_int nrhs);

    /** \brief Codegen bilinear form

        \identifier{su} */
    std::string bilin(const std::string& A, const Sparsity& sp_A,
                      const std::string& x, const std::string& y);

    /** \brief Rank-1 update

        \identifier{sv} */
    std::string rank1(const std::string& A, const Sparsity& sp_A, const std::string& alpha,
                      const std::string& x, const std::string& y);

    /** \brie LogSumExp */
    std::string logsumexp(const std::string& A, casadi_int n);

    /** \brief Multilinear interpolation

        \identifier{sw} */
    std::string interpn(const std::string& res, casadi_int ndim, const std::string& grid,
                        const std::string& offset,
                        const std::string& values, const std::string& x,
                        const std::string& lookup_mode, casadi_int m,
                        const std::string& iw, const std::string& w);

    /** \brief Multilinear interpolation - calculate gradient

        \identifier{sx} */
    std::string interpn_grad(const std::string& grad,
      casadi_int ndim, const std::string& grid,
      const std::string& offset,
      const std::string& values, const std::string& x,
      const std::string& lookup_mode, casadi_int m,
      const std::string& iw, const std::string& w);

    /** \brief Transpose

        \identifier{sy} */
    std::string trans(const std::string& x, const Sparsity& sp_x,
      const std::string& y, const Sparsity& sp_y, const std::string& iw);

    /** \brief QR factorization

        \identifier{sz} */
    std::string qr(const std::string& sp, const std::string& A,
                   const std::string& w, const std::string& sp_v,
                   const std::string& v, const std::string& sp_r,
                   const std::string& r, const std::string& beta,
                   const std::string& prinv, const std::string& pc);

    /** \brief QR solve

        \identifier{t0} */
    std::string qr_solve(const std::string& x, casadi_int nrhs, bool tr,
                         const std::string& sp_v, const std::string& v,
                         const std::string& sp_r, const std::string& r,
                         const std::string& beta, const std::string& prinv,
                         const std::string& pc, const std::string& w);

    /** \\brief LSQR solve

         \identifier{t1} */
    std::string lsqr_solve(const std::string& A, const std::string&x,
                          casadi_int nrhs, bool tr, const std::string& sp, const std::string& w);

    /** \brief LDL factorization

        \identifier{t2} */
    std::string ldl(const std::string& sp_a, const std::string& a,
                   const std::string& sp_lt, const std::string& lt,
                   const std::string& d, const std::string& p,
                   const std::string& w);

    /** \brief LDL solve

        \identifier{t3} */
    std::string ldl_solve(const std::string& x, casadi_int nrhs,
                         const std::string& sp_lt, const std::string& lt,
                         const std::string& d, const std::string& p,
                         const std::string& w);

    /** \brief fmax

        \identifier{t4} */
    std::string fmax(const std::string& x, const std::string& y);

    /** \brief fmin

        \identifier{t5} */
    std::string fmin(const std::string& x, const std::string& y);

    /** \brief mmax

        \identifier{t6} */
    std::string mmax(const std::string& x, casadi_int n, bool is_dense);

    /** \brief mmin

        \identifier{t7} */
    std::string mmin(const std::string& x, casadi_int n, bool is_dense);

    /** \brief vfmax

        \identifier{t8} */
    std::string vfmax(const std::string& x, casadi_int n, const std::string& y);

    /** \brief vfmin

        \identifier{t9} */
    std::string vfmin(const std::string& x, casadi_int n, const std::string& y);

    /** \brief vfmax

        \identifier{ta} */
    std::string vfmax(const std::string& x, const std::string& n, const std::string& y);

    /** \brief vfmin

        \identifier{tb} */
    std::string vfmin(const std::string& x, const std::string& n, const std::string& y);

    /** \brief max

        \identifier{tc} */
    std::string max(const std::string& x, const std::string& y);

    /** \brief min

        \identifier{td} */
    std::string min(const std::string& x, const std::string& y);

    /** \brief norm_inf

        \identifier{te} */
    std::string norm_inf(casadi_int n, const std::string& x);

    /** \brief norm_1

        \identifier{2br} */
    std::string norm_1(casadi_int n, const std::string& x);

    /** 

     * \brief norm_2

                \identifier{256} */
    std::string norm_2(casadi_int n, const std::string& x);

    /** \brief max_viol

        \identifier{tf} */
    std::string max_viol(casadi_int n, const std::string& x,
      const std::string& lb, const std::string& ub);

    /** \brief sum_viol

        \identifier{tg} */
    std::string sum_viol(casadi_int n, const std::string& x,
      const std::string& lb, const std::string& ub);

    /** \brief bound_consistency

        \identifier{th} */
    std::string bound_consistency(casadi_int n, const std::string& x,
      const std::string& lam, const std::string& lbx, const std::string& ubx);

    /** \brief lb_eig

        \identifier{ti} */
    std::string lb_eig(const Sparsity& sp_h, const std::string& h);

    /** \brief regularize

        \identifier{tj} */
    std::string regularize(const Sparsity& sp_h, const std::string& h, const std::string& reg);

    /** \brief convexify

        \identifier{tk} */
    std::string convexify_eval(const ConvexifyData& d,
      const std::string& Hin, const std::string& Hout, const std::string& iw, const std::string& w);

    /** \brief low

        \identifier{tl} */
    std::string low(const std::string& x, const std::string& grid,
      casadi_int ng, casadi_int lookup_mode);

    /** \brief Declare a function

        \identifier{tm} */
    std::string declare(std::string s);

    /** \brief Write a comment line (ignored if not verbose)

        \identifier{tn} */
    void comment(const std::string& s);

    /** \brief Auxiliary functions

        \identifier{to} */
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
      AUX_CLIP_MAX,
      AUX_CLIP_MIN,
      AUX_VECTOR_FMAX,
      AUX_VECTOR_FMIN,
      AUX_NORM_INF,
      AUX_MASKED_NORM_INF,
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
      AUX_QRQP,
      AUX_NLP,
      AUX_SQPMETHOD,
      AUX_FEASIBLESQPMETHOD,
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
      AUX_SUM,
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
      AUX_LOGSUMEXP,
      AUX_SPARSITY,
      AUX_BFGS,
      AUX_ORACLE_CALLBACK,
      AUX_OCP_BLOCK,
      AUX_ORACLE,
      AUX_SCALED_COPY,
      AUX_BLAZING_DE_BOOR,
      AUX_BLAZING_1D_BOOR_EVAL,
      AUX_BLAZING_2D_BOOR_EVAL,
      AUX_BLAZING_3D_BOOR_EVAL,
      AUX_PRINTME,
      AUX_PRINT_SCALAR,
      AUX_PRINT_VECTOR,
      AUX_PRINT_CANONICAL
    };

    /** \brief Add a built-in auxiliary function

        \identifier{tp} */
    void add_auxiliary(Auxiliary f, const std::vector<std::string>& inst = {"casadi_real"});

    /** \brief Add io sparsity patterns of a function

        \identifier{tq} */
    void add_io_sparsities(const std::string& name,
                           const std::vector<Sparsity>& sp_in,
                           const std::vector<Sparsity>& sp_out);

    /** Get work vector name from index */
    std::string work(casadi_int n, casadi_int sz, bool is_ref) const;

    /** Get work vector element from index */
    std::string workel(casadi_int n) const;

    /** \brief Reserve a maximum size of work elements, used for padding of index

        \identifier{2ay} */
    void reserve_work(casadi_int n);

    /** Declare an array */
    static std::string array(const std::string& type, const std::string& name, casadi_int len,
                             const std::string& def=std::string());

    /** \brief  Print casadi_int vector to a c file

        \identifier{tr} */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<casadi_int>& v);

    /** \brief  Print char vector to a c file

        \identifier{27n} */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<char>& v);

    /** \brief  Print real vector to a c file

        \identifier{ts} */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<double>& v);

    /** \brief  Print string vector to a c file

        \identifier{281} */
    void print_vector(std::ostream &s, const std::string& name,
                             const std::vector<std::string>& v);

    /** \brief Print canonical representaion of a matrix

        \identifier{2dk} */
    std::string print_canonical(const Sparsity& sp, const std::string& arg);

    /** \brief Print canonical representaion of a vector

        \identifier{2dl} */
    std::string print_vector(casadi_int sz, const std::string& arg);

    /** \brief Print canonical representaion of a scalar

        \identifier{2dm} */
    std::string print_scalar(const std::string& arg);

    /** \brief Create a copy operation

        \identifier{tt} */
    std::string copy(const std::string& arg, std::size_t n, const std::string& res);
    void copy_check(const std::string& arg, std::size_t n, const std::string& res,
      bool check_lhs=true, bool check_rhs=true);
    void copy_default(const std::string& arg, std::size_t n, const std::string& res,
      const std::string& def,  bool check_rhs=true);

    // Should we elide a copy?
    bool elide_copy(casadi_int sz);

    /** \brief Create a fill operation

        \identifier{tu} */
    std::string fill(const std::string& res, std::size_t n, const std::string& v);

    /** \brief Create a fill operation

        \identifier{tv} */
    std::string clear(const std::string& res, std::size_t n);

    /** \brief Refer to argument

        \identifier{tw} */
    std::string arg(casadi_int i) const;

    /** \brief Refer to resuly

        \identifier{tx} */
    std::string res(casadi_int i) const;

    /** \brief Access thread-local memory

        \identifier{ty} */
    std::string mem(const Function& f);

    /** \brief Sparse assignment

        \identifier{tz} */
    std::string project(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, const Sparsity& sp_res,
                        const std::string& w);

    /** \brief Project triangular part

        \identifier{u0} */
    std::string tri_project(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, bool lower);

    /** \brief Densify

        \identifier{u1} */
    std::string densify(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, bool tr=false);

    /** \brief Sparsify

        \identifier{u2} */
    std::string sparsify(const std::string& arg, const std::string& res,
                         const Sparsity& sp_res, bool tr=false);

    /** \brief Create matrix in MATLAB's MEX format

        \identifier{u3} */
    std::string to_mex(const Sparsity& sp, const std::string& arg);

    /** \brief Get matrix from MATLAB's MEX format

        \identifier{u4} */
    std::string from_mex(std::string& arg,
                         const std::string& res, std::size_t res_off, const Sparsity& sp_res,
                         const std::string& w);

    /** \brief FMU helper functions

        \identifier{257} */
    static std::string fmu_helpers(const std::string& modelname);

    /** \brief Printf

        \identifier{u5} */
    std::string printf(const std::string& str,
                       const std::vector<std::string>& arg=std::vector<std::string>());
    std::string printf(const std::string& str, const std::string& arg1);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2);
    std::string printf(const std::string& str, const std::string& arg1, const std::string& arg2,
                       const std::string& arg3);

    /** \brief Print an operation to a c file

        \identifier{u6} */
    std::string print_op(casadi_int op, const std::string& a0);
    std::string print_op(casadi_int op, const std::string& a0, const std::string& a1);

    /** \brief Slurp a file

        \identifier{u7} */
    std::string file_slurp(const std::string& fname, casadi_int n, const std::string& a);

    /** \brief cache check

        \identifier{u8} */
    std::string cache_check(const std::string& key, const std::string& cache,
        const std::string& loc, casadi_int stride, casadi_int sz, casadi_int key_sz,
        const std::string& val);

    /// Current CasADi version as string
    static std::string casadi_version();

    /// Print file header
    static void stream_open(std::ostream& f, bool cpp);

    /// Print file header
    static void stream_close(std::ostream& f, bool cpp);

    /** \brief Get number of temporary variables needed for all functions

        \identifier{258} */
    void sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const;

  private:

    // Generate casadi_real definition
    void generate_casadi_real(std::ostream &s) const;

    // Generate casadi_int definition
    void generate_casadi_int(std::ostream &s) const;

    // Generate mex entry point
    void generate_mex(std::ostream &s) const;

    // Generate function specific code for Simulink s-Function
    std::string codegen_sfunction(const Function& f) const;

    // Export s-Function to file
    void generate_sfunction(const std::string& name, const std::string& sfunction) const;

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

    // Are we creating a s-function?
    bool with_sfunction;
    std::vector<std::string> added_sfunctions;

    // Unroll arguments?
    bool unroll_args;

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

        \identifier{u9} */
    bool codegen_scalars;

    // Have a flag for exporting/importing symbols
    bool with_export, with_import;

    // Maximum number of declarations per line
    casadi_int max_declarations_per_line;

    // Maximum number of initializer elements per line
    casadi_int max_initializer_elements_per_line;

    // Force the external API to use canonical sparsity
    bool force_canonical;

    // Prefix symbols in DLLs?
    std::string dll_export, dll_import;

    // Prefix
    std::string prefix;

    // std::stringstreams holding the different parts of the file being generated
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

    // Number of zeros/ones
    casadi_int sz_zeros_;
    casadi_int sz_ones_;

    casadi_int padding_length_;

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
    std::multimap<size_t, size_t> added_char_constants_;
    std::multimap<size_t, size_t> added_string_constants_;
    std::map<std::string, std::pair<std::string, std::string> > local_variables_;
    std::map<std::string, std::string> local_default_;
    std::map<const void *, casadi_int> file_scope_double_;
    std::map<const void *, casadi_int> file_scope_integer_;
    std::vector< std::vector<double> > pool_double_defaults_;
    std::map<std::string, casadi_int> pool_double_;

    // Added functions
    struct FunctionMeta {
      // The function object
      Function f;
      // Name in codegen
      std::string codegen_name;
    };
    std::vector<FunctionMeta> added_functions_;

    // Counters for creating unique identifiers
    std::map<std::string, std::map<FunctionInternal*, casadi_int> > added_wrappers_;

    // Constants
    std::vector<std::vector<double> > double_constants_;
    std::vector<std::vector<casadi_int> > integer_constants_;
    std::vector<std::vector<char> > char_constants_;
    std::vector<std::vector<std::string> > string_constants_;

    // Does any function need thread-local memory?
    bool needs_mem_;

    // Hash a vector
    static size_t hash(const std::vector<double>& v);
    static size_t hash(const std::vector<casadi_int>& v);
    static size_t hash(const std::vector<char>& v);
    static size_t hash(const std::vector<std::string>& v);

    std::string wrapper(const Function& base, const std::string& name);

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
