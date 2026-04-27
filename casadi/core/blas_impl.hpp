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

namespace casadi { class CodeGenerator; }


/// \cond INTERNAL
namespace casadi {

  // CBLAS_TRANSPOSE values, mirrored to avoid pulling in <cblas.h>.
  // Plugins built against OpenBLAS/BLASFEO can pass these straight through.
  enum {
    CASADI_BLAS_NO_TRANS = 111,
    CASADI_BLAS_TRANS    = 112
  };

  /** \brief Blas interface

      Pluggable dense matrix-multiply backend. Layout is column-major.
      The dispatch entry Blas::dgemm(shorthand, ...) is keyed by a small
      integer assigned at plugin registration time; shorthand 0 is
      reserved for the built-in "reference" plugin and takes a fast path
      that bypasses the function-pointer indirection.

      Plugins registered via the standard PluginInterface mechanism
      expose a single dgemm function pointer in Plugin::Exposed. External
      plugins (openblas, blasfeo) live in their own shared libraries and
      are dlopen'd on first use.

      \author Joris Gillis
      \date 2026

      */
  class CASADI_EXPORT
  Blas : public PluginInterface<Blas> {
  public:
    /// CBLAS-style dgemm: C := alpha*op(A)*op(B) + beta*C
    /// Layout is column-major; transa/transb are CASADI_BLAS_NO_TRANS / CASADI_BLAS_TRANS.
    typedef void (* Dgemm)(int transa, int transb,
                           casadi_int m, casadi_int n, casadi_int k,
                           double alpha,
                           const double* A, casadi_int lda,
                           const double* B, casadi_int ldb,
                           double beta,
                           double* C, casadi_int ldc);

    /// Codegen counterpart of `mtimes`: emit a C statement that performs
    /// C += A*B (column-major, NoTrans, contiguous). The plugin owns its
    /// own emission — extern decls, helper #ifdefs, dimension layout —
    /// so plugin-specific knowledge stays out of libcasadi core.
    typedef void (* CodegenMtimes)(CodeGenerator& g,
                                   const std::string& A,
                                   casadi_int m, casadi_int k,
                                   const std::string& B, casadi_int n,
                                   const std::string& C);

    // Creator function (unused; kept for PluginInterface uniformity).
    typedef Blas* (*Creator)();

    static const std::string meta_doc;

    struct Exposed {
      Dgemm dgemm;
      CodegenMtimes codegen_mtimes;
    };

    /// Collection of registered plugins, keyed by name.
    static std::map<std::string, Plugin> solvers_;

    /// O(1) dispatch by shorthand. Index 0 is reserved for "reference".
    /// Pre-reserved to MAX_PLUGINS so push_back never reallocates;
    /// element pointers therefore stay stable for process lifetime.
    static std::vector<const Plugin*> dispatch_;

    /// Hard cap on the number of BLAS plugins per process.
    static constexpr unsigned char MAX_PLUGINS = 255;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix used by PluginInterface to derive DLL / register-symbol names.
    static const std::string infix_;

    /// Resolve a plugin name to its dispatch shorthand.
    /// "reference" is hard-coded to shorthand 0 (no Plugin entry, no DLL).
    /// External plugin names auto-load on first use.
    /// Intended to be called once per Multiplication node, at construction.
    static unsigned char shorthand_for(const std::string& name);

    /// Full CBLAS-shape dispatch: C := alpha*op(A)*op(B) + beta*C, column-major.
    /// For shorthand==0 (reference), inlines a libcasadi-internal impl that
    /// fast-paths the canonical (alpha=1, beta=1, contiguous) case to
    /// casadi_mtimes_dense. Other shorthands trigger an indirect call
    /// through dispatch_[shorthand]->exposed.dgemm.
    static void dgemm(unsigned char shorthand,
                      int transa, int transb,
                      casadi_int m, casadi_int n, casadi_int k,
                      double alpha,
                      const double* A, casadi_int lda,
                      const double* B, casadi_int ldb,
                      double beta,
                      double* C, casadi_int ldc);

    /// Canonical dense matrix-multiply-accumulate: C += A * B,
    /// column-major, both factors NoTrans, contiguous (lda=m, ldb=k, ldc=m).
    /// This is the runtime counterpart of codegen_mtimes(): same call shape
    /// inside Multiplication eval and generate.
    static void mtimes(unsigned char shorthand,
                       const double* A, casadi_int m, casadi_int k,
                       const double* B, casadi_int n,
                       double* C);

    /// Built-in dense dgemm (used internally by the shorthand==0 fast path).
    /// Reachable directly for unit tests; normal callers should use mtimes/dgemm.
    static void reference_dgemm(int transa, int transb,
                             casadi_int m, casadi_int n, casadi_int k,
                             double alpha,
                             const double* A, casadi_int lda,
                             const double* B, casadi_int ldb,
                             double beta,
                             double* C, casadi_int ldc);

    /// Reverse lookup: shorthand -> plugin name. Returns "reference" for 0;
    /// otherwise the registered plugin's name (caller must ensure `shorthand`
    /// was previously returned by shorthand_for()).
    static const char* name_for_shorthand(unsigned char shorthand);

    /// Emit a C code line that performs C += A*B (column-major dense),
    /// dispatched to the requested BLAS plugin. Shorthand 0 (reference)
    /// emits the built-in casadi_mtimes_dense call; non-reference branches
    /// will be wired up alongside their codegen support. Codegen counterpart
    /// of Blas::mtimes().
    static void codegen_mtimes(CodeGenerator& g, unsigned char shorthand,
                               const std::string& A,
                               casadi_int m, casadi_int k,
                               const std::string& B, casadi_int n,
                               const std::string& C);

  };

} // namespace casadi

/// \endcond

#endif // CASADI_BLAS_IMPL_HPP
