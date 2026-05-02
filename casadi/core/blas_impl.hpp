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
    typedef void   (* Dcopy)(const double* x, casadi_int n, double* y);

    /* \brief Codegen counterpart of an L1 op: emit the auxiliary block */
    typedef void (* CodegenL1Aux)(CodeGenerator& g,
                                  const std::vector<std::string>& inst);

    /* \brief Codegen counterpart of an L3 op (mtimes_dense): emit the
     * auxiliary block.  Same shape as the L1 aux hook so a BLAS plugin
     * can override the canonical reference loop with one that calls
     * cblas_dgemm (or any equivalent). */
    typedef void (* CodegenL3Aux)(CodeGenerator& g,
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
      Dcopy dcopy;
      CodegenL1Aux codegen_axpy_aux;
      CodegenL1Aux codegen_dot_aux;
      CodegenL1Aux codegen_scal_aux;
      CodegenL1Aux codegen_nrm2_aux;
      CodegenL1Aux codegen_asum_aux;
      /* \brief Optional override for AUX_MTIMES_DENSE: when set, the
       * plugin emits a custom casadi_mtimes_dense body (e.g. cblas_dgemm
       * wrapper).  When null, codegen falls back to the reference loop. */
      CodegenL3Aux codegen_mtimes_dense_aux;
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

    /* \brief Canonical contiguous mtimes-accumulate: C += op(A)*B
     *  tr=0: A is m×k, B is k×n, C is m×n.
     *  tr=1: op(A)=A^T; A is m×k, B is m×n, C is k×n. */
    static void mtimes(casadi_int shorthand,
                       const double* A, casadi_int m, casadi_int k,
                       const double* B, casadi_int n,
                       double* C, casadi_int tr = 0);


    /* \brief Emit a C statement that performs C += A*B; codegen counterpart of mtimes */
    static void codegen_mtimes(CodeGenerator& g, casadi_int shorthand,
                               const std::string& A,
                               casadi_int m, casadi_int k,
                               const std::string& B, casadi_int n,
                               const std::string& C);

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
    /* \brief Try-emit aux block for casadi_mtimes_dense via active plugin
     * (false on fallback to the reference loop). */
    static bool codegen_mtimes_dense_aux(CodeGenerator& g,
                                         const std::vector<std::string>& inst);
  };

} // namespace casadi

/// \endcond

#endif // CASADI_BLAS_IMPL_HPP
