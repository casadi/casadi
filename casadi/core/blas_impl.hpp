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


#ifndef CASADI_BLAS_IMPL_HPP
#define CASADI_BLAS_IMPL_HPP

#include "blas.hpp"
#include "plugin_interface.hpp"
#include "runtime/casadi_runtime.hpp"

#include <cstring>

namespace casadi { class CodeGenerator; }


/// \cond INTERNAL
namespace casadi {

  /* \brief CBLAS_TRANSPOSE values, mirrored so plugins don't need <cblas.h> */
  enum {
    CASADI_BLAS_NO_TRANS = 111,
    CASADI_BLAS_TRANS    = 112
  };

  /** \brief Pluggable dense BLAS backend, dispatched per-call by shorthand

      Shorthand 0 is the built-in reference impl (fast path, no indirection).
      External plugins (classic, blasfeo, ...) live in their own DLLs and are
      dlopen'd on first use; each populates Exposed with its own function
      pointers.

      \author Joris Gillis
      \date 2026

      \identifier{2ge} */
  class CASADI_EXPORT
  Blas : public PluginInterface<Blas> {
  public:
    /* \brief CBLAS-style dgemm: C := alpha*op(A)*op(B) + beta*C; column-major */
    typedef void (* Dgemm)(int transa, int transb,
                           casadi_int m, casadi_int n, casadi_int k,
                           double alpha,
                           const double* A, casadi_int lda,
                           const double* B, casadi_int ldb,
                           double beta,
                           double* C, casadi_int ldc);

    /* \brief Codegen counterpart of mtimes: emit C += A*B at the call site */
    typedef void (* CodegenMtimes)(CodeGenerator& g,
                                   const std::string& A,
                                   casadi_int m, casadi_int k,
                                   const std::string& B, casadi_int n,
                                   const std::string& C);

    /* \brief L1 runtime hook signatures (double-only; templates handle other T) */
    typedef void   (* Daxpy)(casadi_int n, double alpha, const double* x, double* y);
    typedef double (* Ddot )(casadi_int n, const double* x, const double* y);
    typedef void   (* Dscal)(casadi_int n, double alpha, double* x);
    typedef double (* Dnrm2)(casadi_int n, const double* x);
    typedef double (* Dasum)(casadi_int n, const double* x);

    /* \brief Codegen counterpart of an L1 op: emit the auxiliary block */
    typedef void (* CodegenL1Aux)(CodeGenerator& g,
                                  const std::vector<std::string>& inst);

    /* \brief Creator function (unused; kept for PluginInterface uniformity) */
    typedef Blas* (*Creator)();

    static const std::string meta_doc;

    /* \brief Function-pointer table populated by each registered plugin */
    struct Exposed {
      Dgemm dgemm;
      CodegenMtimes codegen_mtimes;
      Daxpy daxpy;
      Ddot  ddot;
      Dscal dscal;
      Dnrm2 dnrm2;
      Dasum dasum;
      CodegenL1Aux codegen_axpy_aux;
      CodegenL1Aux codegen_dot_aux;
      CodegenL1Aux codegen_scal_aux;
      CodegenL1Aux codegen_nrm2_aux;
      CodegenL1Aux codegen_asum_aux;
    };

    /* \brief Collection of registered plugins, keyed by name */
    static std::map<std::string, Plugin> solvers_;

    /* \brief Fast lookup structure for BLAS plugins */
    static std::vector<const Plugin*> dispatch_;

    /* \brief Active default BLAS shorthand (index into dispatch_; 0 = reference) */
    static casadi_int default_;

    /* \brief Get the active default BLAS name */
    static std::string getDefault();

    /* \brief Resolve a plugin name to its dispatch shorthand */
    static casadi_int shorthand_for(const std::string& name);

    /* \brief Reverse lookup: shorthand -> plugin name ("reference" for 0) */
    static const char* name_for_shorthand(casadi_int shorthand);
  private:
    friend class GlobalOptions;

    /* \brief Set default_; private, callers go through GlobalOptions::setDefaultBlas */
    static void setDefault(const std::string& name);

  public:

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /* \brief Infix used by PluginInterface to derive DLL / register-symbol names */
    static const std::string infix_;

    /* \brief Full CBLAS dgemm dispatched by shorthand (column-major) */
    static void dgemm(casadi_int shorthand,
                      int transa, int transb,
                      casadi_int m, casadi_int n, casadi_int k,
                      double alpha,
                      const double* A, casadi_int lda,
                      const double* B, casadi_int ldb,
                      double beta,
                      double* C, casadi_int ldc);

    /* \brief Built-in dense dgemm (shorthand-0 fast path; exposed for tests) */
    static void reference_dgemm(int transa, int transb,
                             casadi_int m, casadi_int n, casadi_int k,
                             double alpha,
                             const double* A, casadi_int lda,
                             const double* B, casadi_int ldb,
                             double beta,
                             double* C, casadi_int ldc);

    /* \brief Canonical contiguous mtimes-accumulate: C += A*B */
    static void mtimes(casadi_int shorthand,
                       const double* A, casadi_int m, casadi_int k,
                       const double* B, casadi_int n,
                       double* C);

    /* \brief Emit a C statement that performs C += A*B; codegen counterpart of mtimes */
    static void codegen_mtimes(CodeGenerator& g, casadi_int shorthand,
                               const std::string& A,
                               casadi_int m, casadi_int k,
                               const std::string& B, casadi_int n,
                               const std::string& C);

    /* \brief y <- x (zero-fill if x==NULL) */
    template<typename T>
    static void copy(const T* x, casadi_int n, T* y) { casadi_copy(x, n, y); }

    /* \brief y += alpha*x */
    template<typename T>
    static void axpy(casadi_int n, T alpha, const T* x, T* y) {
      casadi_axpy(n, alpha, x, y);
    }

    /* \brief dot(x, y) */
    template<typename T>
    static T dot(casadi_int n, const T* x, const T* y) {
      return casadi_dot(n, x, y);
    }

    /* \brief x *= alpha */
    template<typename T>
    static void scal(casadi_int n, T alpha, T* x) {
      casadi_scal(n, alpha, x);
    }

    /* \brief 2-norm of x */
    template<typename T>
    static T norm_2(casadi_int n, const T* x) {
      return casadi_norm_2(n, x);
    }

    /* \brief 1-norm of x */
    template<typename T>
    static T norm_1(casadi_int n, const T* x) {
      return casadi_norm_1(n, x);
    }

    /* \brief Emit aux block for casadi_copy (memcpy/memset; never plugin-dispatched) */
    static void codegen_copy_aux(CodeGenerator& g,
                                 const std::vector<std::string>& inst);
    /* \brief Try-emit aux block for casadi_axpy via active plugin (false on fallback) */
    static bool codegen_axpy_aux(CodeGenerator& g,
                                 const std::vector<std::string>& inst);
    /* \brief Try-emit aux block for casadi_dot via active plugin (false on fallback) */
    static bool codegen_dot_aux(CodeGenerator& g,
                                const std::vector<std::string>& inst);
    /* \brief Try-emit aux block for casadi_scal via active plugin (false on fallback) */
    static bool codegen_scal_aux(CodeGenerator& g,
                                 const std::vector<std::string>& inst);
    /* \brief Try-emit aux block for casadi_norm_2 via active plugin (false on fallback) */
    static bool codegen_norm_2_aux(CodeGenerator& g,
                                   const std::vector<std::string>& inst);
    /* \brief Try-emit aux block for casadi_norm_1 via active plugin (false on fallback) */
    static bool codegen_norm_1_aux(CodeGenerator& g,
                                   const std::vector<std::string>& inst);
  };

  template<>
  inline void Blas::copy<double>(const double* x, casadi_int n, double* y) {
    if (!y) return;
    if (x) std::memcpy(y, x, n * sizeof(double));
    else   std::memset(y, 0, n * sizeof(double));
  }

  template<>
  inline void Blas::axpy<double>(casadi_int n, double alpha,
                                 const double* x, double* y) {
    if (!default_) return casadi_axpy(n, alpha, x, y);
    if (!x || !y) return;
    Daxpy fn = dispatch_[default_]->exposed.daxpy;
    if (fn) return fn(n, alpha, x, y);
    return casadi_axpy(n, alpha, x, y);
  }

  template<>
  inline double Blas::dot<double>(casadi_int n, const double* x, const double* y) {
    if (!default_) return casadi_dot(n, x, y);
    Ddot fn = dispatch_[default_]->exposed.ddot;
    if (fn) return fn(n, x, y);
    return casadi_dot(n, x, y);
  }

  template<>
  inline void Blas::scal<double>(casadi_int n, double alpha, double* x) {
    if (!default_) return casadi_scal(n, alpha, x);
    if (!x) return;
    Dscal fn = dispatch_[default_]->exposed.dscal;
    if (fn) return fn(n, alpha, x);
    return casadi_scal(n, alpha, x);
  }

  template<>
  inline double Blas::norm_2<double>(casadi_int n, const double* x) {
    if (!default_) return casadi_norm_2(n, x);
    Dnrm2 fn = dispatch_[default_]->exposed.dnrm2;
    if (fn) return fn(n, x);
    return casadi_norm_2(n, x);
  }

  template<>
  inline double Blas::norm_1<double>(casadi_int n, const double* x) {
    if (!default_) return casadi_norm_1(n, x);
    if (!x) return 0;
    Dasum fn = dispatch_[default_]->exposed.dasum;
    if (fn) return fn(n, x);
    return casadi_norm_1(n, x);
  }

} // namespace casadi

/// \endcond

#endif // CASADI_BLAS_IMPL_HPP
