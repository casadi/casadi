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


#include "blas_impl.hpp"
#include "code_generator.hpp"
#include "global_options.hpp"
#include "runtime/casadi_runtime.hpp"

#include <cstring>

namespace casadi {

std::map<std::string, Blas::Plugin> Blas::solvers_;
std::vector<const Blas::Plugin*> Blas::dispatch_;

// Default-plugin storage is split:
//   - Blas::default_              (casadi_int): hot-path read by every L1
//                                 dispatcher; integer shorthand. Defined here.
//   - GlobalOptions::default_blas_ (std::string): user-facing canonical name.
//                                  Defined in global_options.cpp.
// GlobalOptions::setDefaultBlas keeps both in sync; Blas::setDefault is the
// internal entry point that updates only the integer shorthand.
casadi_int Blas::default_ = 0;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
std::mutex Blas::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

const std::string Blas::infix_ = "blas";

void Blas::reference_dgemm(int transa, int transb,
                        casadi_int m, casadi_int n, casadi_int k,
                        double alpha,
                        const double* A, casadi_int lda,
                        const double* B, casadi_int ldb,
                        double beta,
                        double* C, casadi_int ldc) {
  const bool transa_yes = (transa == CASADI_BLAS_TRANS);
  const bool transb_yes = (transb == CASADI_BLAS_TRANS);

  // Fast path: canonical contiguous column-major layout, alpha=beta=1, transb=No.
  // casadi_mtimes_dense takes (x, nrow_x, ncol_x, y, ncol_y, z, tr) and
  // computes z += op(x) * y where:
  //   tr=0: x is (nrow_x, ncol_x), y is (ncol_x, ncol_y) -> z (nrow_x, ncol_y)
  //   tr=1: x is (nrow_x, ncol_x), y is (nrow_x, ncol_y) -> z (ncol_x, ncol_y)
  const casadi_int lda_canonical = transa_yes ? k : m;
  const bool canonical = !transb_yes
      && alpha == 1.0 && beta == 1.0
      && lda == lda_canonical && ldb == k && ldc == m;

  if (canonical) {
    if (transa_yes) {
      // z (m, n) += A^T (m, k from k-by-m) * B (k, n)
      casadi_mtimes_dense<double>(A, k, m, B, n, C, 1);
    } else {
      // z (m, n) += A (m, k) * B (k, n)
      casadi_mtimes_dense<double>(A, m, k, B, n, C, 0);
    }
    return;
  }

  // General fallback: scale C, then accumulate.
  if (beta == 0.0) {
    for (casadi_int j = 0; j < n; ++j) {
      double* col = C + j * ldc;
      for (casadi_int i = 0; i < m; ++i) col[i] = 0.0;
    }
  } else if (beta != 1.0) {
    for (casadi_int j = 0; j < n; ++j) {
      double* col = C + j * ldc;
      for (casadi_int i = 0; i < m; ++i) col[i] *= beta;
    }
  }

  for (casadi_int j = 0; j < n; ++j) {
    for (casadi_int l = 0; l < k; ++l) {
      const double b_lj = transb_yes ? B[j + l * ldb] : B[l + j * ldb];
      if (b_lj == 0.0) continue;
      const double scl = alpha * b_lj;
      double* col = C + j * ldc;
      if (transa_yes) {
        const double* a_col = A + l;  // walks A^T's column l = A's row l
        // a_il = A[l + i * lda]
        for (casadi_int i = 0; i < m; ++i) col[i] += scl * a_col[i * lda];
      } else {
        const double* a_col = A + l * lda;  // A's column l
        for (casadi_int i = 0; i < m; ++i) col[i] += scl * a_col[i];
      }
    }
  }
}

// Shorthand 0 is reserved for the built-in reference impl. It has no Plugin
// entry; dispatch_[0] is a permanent nullptr sentinel that's never read
// because every dispatch path special-cases sh==0.
static const char* REFERENCE_DOC =
    "Built-in dense BLAS implementation, no external dependency. "
    "Used by default and as the fallback when other plugins are unavailable.";

casadi_int Blas::shorthand_for(const std::string& name) {
  if (name == "reference") return 0;
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::lock_guard<std::mutex> lock(Blas::mutex_solvers_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

  // Lazy-init the dispatch_ vector: plant the slot-0 sentinel.
  if (dispatch_.empty()) {
    dispatch_.push_back(nullptr);  // index 0 == reference, never indexed
  }

  auto it = solvers_.find(name);
  if (it == solvers_.end()) {
    // Auto-load (lock already held, so pass needs_lock=false)
    load_plugin(name, true, false);
    it = solvers_.find(name);
    casadi_assert_dev(it != solvers_.end());
  }

  // Find existing shorthand by pointer identity (small N, linear scan is fine).
  // Start at 1 — slot 0 is the reference sentinel.
  const Plugin* p = &it->second;
  for (casadi_int sh = 1; sh < static_cast<casadi_int>(dispatch_.size()); ++sh) {
    if (dispatch_[sh] == p) return sh;
  }

  // Newly registered plugin: assign next shorthand
  dispatch_.push_back(p);
  return static_cast<casadi_int>(dispatch_.size() - 1);
}

void Blas::dgemm(casadi_int shorthand,
                 int transa, int transb,
                 casadi_int m, casadi_int n, casadi_int k,
                 double alpha,
                 const double* A, casadi_int lda,
                 const double* B, casadi_int ldb,
                 double beta,
                 double* C, casadi_int ldc) {
  if (shorthand == 0) {
    // Reference fast path: skip indirection entirely.
    reference_dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    return;
  }
  // External plugin: trust the caller obtained `shorthand` via shorthand_for(),
  // which guarantees dispatch_[shorthand] is populated. Hot path, no lock.
  casadi_assert_dev(shorthand < static_cast<casadi_int>(dispatch_.size()));
  casadi_assert_dev(dispatch_[shorthand] != nullptr);
  dispatch_[shorthand]->exposed.dgemm(
      transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

void Blas::mtimes(casadi_int shorthand,
                  const double* A, casadi_int m, casadi_int k,
                  const double* B, casadi_int n,
                  double* C) {
  if (shorthand == 0) {
    // Reference fast path: arguments are already canonical here, so skip
    // dgemm's canonical-detection branch and call the loop directly.
    casadi_mtimes_dense<double>(A, m, k, B, n, C, /*tr=*/0);
    return;
  }
  // External plugin: indirect through dispatch_, no lock. The caller
  // obtained `shorthand` via shorthand_for() so the slot is populated.
  casadi_assert_dev(shorthand < static_cast<casadi_int>(dispatch_.size()));
  casadi_assert_dev(dispatch_[shorthand] != nullptr);
  dispatch_[shorthand]->exposed.dgemm(
      CASADI_BLAS_NO_TRANS, CASADI_BLAS_NO_TRANS,
      m, n, k, 1.0, A, m, B, k, 1.0, C, m);
}

const char* Blas::name_for_shorthand(casadi_int shorthand) {
  if (shorthand == 0) return "reference";
  // Read-only fast path: dispatch_ entries are stable for process lifetime
  // and shorthand_for() is the sole writer. No lock needed; trust caller.
  casadi_assert_dev(shorthand < static_cast<casadi_int>(dispatch_.size()));
  casadi_assert_dev(dispatch_[shorthand] != nullptr);
  return dispatch_[shorthand]->name;
}

void Blas::codegen_mtimes(CodeGenerator& g, casadi_int shorthand,
                          const std::string& A,
                          casadi_int m, casadi_int k,
                          const std::string& B, casadi_int n,
                          const std::string& C) {
  if (shorthand == 0) {
    // Reference fast path: emit the built-in casadi_mtimes_dense call.
    g << g.mtimes(A, m, k, B, n, C, false) << '\n';
    return;
  }
  // External plugin: trust the caller obtained `shorthand` via shorthand_for(),
  // so dispatch_[shorthand] is populated and its codegen_mtimes is non-null
  // (otherwise the plugin would have failed to register).
  casadi_assert_dev(shorthand < static_cast<casadi_int>(dispatch_.size()));
  casadi_assert_dev(dispatch_[shorthand] != nullptr);
  const Plugin* p = dispatch_[shorthand];
  casadi_assert(p->exposed.codegen_mtimes != nullptr,
      "BLAS plugin '" + std::string(p->name) + "' does not implement codegen.");
  p->exposed.codegen_mtimes(g, A, m, k, B, n, C);
}

bool has_blas(const std::string& name) {
  if (name == "reference") return true;
  return Blas::has_plugin(name);
}

void load_blas(const std::string& name) {
  if (name == "reference") return;
  Blas::load_plugin(name);
}

std::string doc_blas(const std::string& name) {
  if (name == "reference") return REFERENCE_DOC;
  return Blas::getPlugin(name).doc;
}

// Internal entry point: sets the integer shorthand only. shorthand_for()
// validates the name and lazy-loads the plugin DLL if needed. The matching
// GlobalOptions::setDefaultBlas (in global_options.cpp) is the public entry
// point that mirrors the name string and calls this.
void Blas::setDefault(const std::string& name) {
  default_ = shorthand_for(name);
}

std::string Blas::getDefault() {
  return name_for_shorthand(default_);
}

void Blas::codegen_copy_aux(CodeGenerator& g,
                            const std::vector<std::string>& inst) {
  g.add_include("string.h");
  g.auxiliaries << g.sanitize_source(
      "// SYMBOL \"copy\"\n"
      "void casadi_copy(const casadi_real* x, casadi_int n, casadi_real* y) {\n"
      "  if (!y) return;\n"
      "  if (x) memcpy(y, x, n*sizeof(casadi_real));\n"
      "  else   memset(y, 0, n*sizeof(casadi_real));\n"
      "}\n",
      inst);
}

bool Blas::codegen_axpy_aux(CodeGenerator& g,
                            const std::vector<std::string>& inst) {
  if (!default_) return false;
  CodegenL1Aux fn = dispatch_[default_]->exposed.codegen_axpy_aux;
  if (!fn) return false;
  fn(g, inst);
  return true;
}

bool Blas::codegen_dot_aux(CodeGenerator& g,
                           const std::vector<std::string>& inst) {
  if (!default_) return false;
  CodegenL1Aux fn = dispatch_[default_]->exposed.codegen_dot_aux;
  if (!fn) return false;
  fn(g, inst);
  return true;
}

bool Blas::codegen_scal_aux(CodeGenerator& g,
                            const std::vector<std::string>& inst) {
  if (!default_) return false;
  CodegenL1Aux fn = dispatch_[default_]->exposed.codegen_scal_aux;
  if (!fn) return false;
  fn(g, inst);
  return true;
}

bool Blas::codegen_norm_2_aux(CodeGenerator& g,
                              const std::vector<std::string>& inst) {
  if (!default_) return false;
  CodegenL1Aux fn = dispatch_[default_]->exposed.codegen_nrm2_aux;
  if (!fn) return false;
  fn(g, inst);
  return true;
}

bool Blas::codegen_norm_1_aux(CodeGenerator& g,
                              const std::vector<std::string>& inst) {
  if (!default_) return false;
  CodegenL1Aux fn = dispatch_[default_]->exposed.codegen_asum_aux;
  if (!fn) return false;
  fn(g, inst);
  return true;
}

} // namespace casadi
