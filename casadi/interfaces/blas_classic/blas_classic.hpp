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


#ifndef CASADI_BLAS_CLASSIC_HPP
#define CASADI_BLAS_CLASSIC_HPP

#include "casadi/core/blas_impl.hpp"
#include <casadi/interfaces/blas_classic/casadi_blas_classic_export.h>

extern "C" {
  /// Fortran BLAS dgemm (column-major, 32-bit ints).
  /// Available from any BLAS library configured into the CasADi build
  /// via the `lapack` CMake target (OpenBLAS, MKL, system BLAS, ...).
  void dgemm_(const char* transa, const char* transb,
              const int* m, const int* n, const int* k,
              const double* alpha,
              const double* a, const int* lda,
              const double* b, const int* ldb,
              const double* beta,
              double* c, const int* ldc);

  /// Fortran BLAS L1 (column-major, 32-bit ints, scalars passed by reference).
  void daxpy_(const int* n, const double* alpha,
              const double* x, const int* incx,
              double* y, const int* incy);
  double ddot_(const int* n,
               const double* x, const int* incx,
               const double* y, const int* incy);
  void dscal_(const int* n, const double* alpha,
              double* x, const int* incx);
  double dnrm2_(const int* n, const double* x, const int* incx);
  double dasum_(const int* n, const double* x, const int* incx);
}

/** \defgroup plugin_Blas_classic Title
    \par
*/

/** \pluginsection{Blas,classic} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{Blas,classic}

    Wraps the BLAS Fortran ABI (`dgemm_`) of whatever library was
    configured into the CasADi build (typically OpenBLAS, but any
    Fortran-ABI BLAS will do).

    @copydoc Blas_doc
    @copydoc plugin_Blas_classic
    \author Joris Gillis
    \date 2026
  */
  class CASADI_BLAS_CLASSIC_EXPORT BlasClassic : public Blas {
  public:
    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BLAS_CLASSIC_HPP
