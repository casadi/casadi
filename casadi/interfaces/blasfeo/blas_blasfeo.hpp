/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_BLAS_BLASFEO_HPP
#define CASADI_BLAS_BLASFEO_HPP

#include "casadi/core/blas_impl.hpp"
#include <casadi/interfaces/blasfeo/casadi_blas_blasfeo_export.h>

extern "C" {
  /// BLASFEO's namespaced Fortran-ABI dgemm.
  /// Non-const pointers to match the BLASFEO header exactly.
  void blasfeo_blas_dgemm(char* ta, char* tb,
                          int* m, int* n, int* k,
                          double* alpha,
                          double* A, int* lda,
                          double* B, int* ldb,
                          double* beta,
                          double* C, int* ldc);
}

/** \defgroup plugin_Blas_blasfeo Title
    \par

    \identifier{2g9} */

/** \pluginsection{Blas,blasfeo} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{Blas,blasfeo}

    Wraps BLASFEO's `blasfeo_blas_dgemm` (namespaced Fortran ABI).
    BLASFEO is tuned for small/medium dense matrices common in
    embedded optimization workloads.

    @copydoc Blas_doc
    @copydoc plugin_Blas_blasfeo
    \author Joris Gillis
    \date 2026
  */
  class CASADI_BLAS_BLASFEO_EXPORT BlasBlasfeo : public Blas {
  public:
    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BLAS_BLASFEO_HPP
