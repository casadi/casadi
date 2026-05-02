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


#ifndef CASADI_BLAS_HPP
#define CASADI_BLAS_HPP

#include <casadi/core/casadi_export.h>
#include <casadi/core/casadi_types.hpp>
#include <string>

namespace casadi {

  /** \defgroup main_blas Title
      \par

      Pluggable dense matrix-multiply backend. The built-in "reference"
      plugin requires no external dependency; "classic" and "blasfeo"
      plugins (when built) load lazily from their respective shared
      libraries on first use.

      \generalsection{Blas}
      \pluginssection{Blas}

      \author Joris Gillis
      \date 2026

      \identifier{2gb} */

  /* \brief Check if a particular BLAS plugin is available */
  CASADI_EXPORT bool has_blas(const std::string& name);

  /* \brief Explicitly load a BLAS plugin dynamically */
  CASADI_EXPORT void load_blas(const std::string& name);

  /* \brief Get the documentation string for a BLAS plugin */
  CASADI_EXPORT std::string doc_blas(const std::string& name);

#ifndef SWIG
  /* \brief Vm-side BLAS dispatch wrapper around casadi_mtimes_dense.
   *
   *  Distinct symbol from casadi_mtimes_dense: the <double>
   *  specialization (defined out-of-line in blas.cpp) routes through
   *  Blas::mtimes(default_, ...) -- the integration point that makes
   *  setDefaultBlas("classic"|"blasfeo"|...) flip vm-side dispatch in
   *  callers like casadi_condensing_eval<double>.
   *
   *  Why a separate name rather than specializing casadi_mtimes_dense
   *  directly?  The spec body must call the unspecialized reference
   *  primary on the shorthand-0 fast path (Blas::mtimes(0, ...) calls
   *  casadi_mtimes_dense<double> -- the triple loop).  Sharing the
   *  symbol name would recurse.
   *
   *  Visibility contract: every TU that instantiates a templated
   *  function whose body calls casadi_blas_mtimes (currently
   *  only casadi_condensing_eval<double>, instantiated in conic.cpp)
   *  MUST include this header BEFORE any chain that parses the
   *  template body.  Phase-1 lookup at template-parse time has to see
   *  the declaration; otherwise phase-2 ADL on primitive pointer
   *  arguments finds nothing and the compile errors out at the call
   *  site.  Forgetting the include is loud, never silent.
   *
   *  In codegen output, casadi_condensing.hpp's call sites are
   *  rewritten via C-REPLACE to plain casadi_mtimes_dense (which the
   *  plugin AUX hook overrides), so the wrapper has no codegen
   *  footprint. */
  template<typename T1>
  void casadi_blas_mtimes(const T1* x, casadi_int nrow_x, casadi_int ncol_x,
      const T1* y, casadi_int ncol_y, T1* z, casadi_int tr) {
    casadi_mtimes_dense(x, nrow_x, ncol_x, y, ncol_y, z, tr);
  }

  template<>
  void CASADI_EXPORT casadi_blas_mtimes<double>(
      const double* A, casadi_int m, casadi_int k,
      const double* B, casadi_int n, double* C, casadi_int tr);
#endif // SWIG

} // namespace casadi


#endif // CASADI_BLAS_HPP
