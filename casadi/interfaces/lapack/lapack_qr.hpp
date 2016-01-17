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


#ifndef CASADI_LAPACK_QR_HPP
#define CASADI_LAPACK_QR_HPP

#include "casadi/core/function/linsol.hpp"
#include <casadi/interfaces/lapack/casadi_linsol_lapackqr_export.h>

/** \defgroup plugin_Linsol_lapackqr
*
* This class solves the linear system <tt>A.x=b</tt> by making an QR factorization of A: \n
* <tt>A = Q.R</tt>, with Q orthogonal and R upper triangular
*/

/** \pluginsection{Linsol,lapackqr} */

/// \cond INTERNAL
namespace casadi {

  /// QR-factorize dense matrix (lapack)
  extern "C" void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
                          double *work, int *lwork, int *info);

  /// Multiply right hand side with Q-transpose (lapack)
  extern "C" void dormqr_(char *side, char *trans, int *n, int *m, int *k, double *a,
                          int *lda, double *tau, double *c, int *ldc,
                          double *work, int *lwork, int *info);

  /// Solve upper triangular system (lapack)
  extern "C" void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
                         double *alpha, double *a, int *lda, double *b, int *ldb);

  /** \brief  \pluginbrief{Linsol,lapackqr}
   *
   @copydoc Linsol_doc
   @copydoc plugin_Linsol_lapackqr
   *
   */
  class CASADI_LINSOL_LAPACKQR_EXPORT LapackQr : public Linsol {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LapackQr(const std::string& name, const Sparsity& sparsity, int nrhs);

    /** \brief  Create a new Linsol */
    static Linsol* creator(const std::string& name, const Sparsity& sp, int nrhs) {
      return new LapackQr(name, sp, nrhs);
    }

    // Destructor
    virtual ~LapackQr();

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Allocate memory block */
    virtual Memory* memory() const;

    // Factorize the linear system
    virtual void linsol_factorize(Memory& mem, const double* A) const;

    // Solve the linear system
    virtual void linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    virtual const char* plugin_name() const { return "lapackqr";}
  };

  struct CASADI_LINSOL_LAPACKQR_EXPORT LapackQrMemory : public Memory {
    // Destructor
    virtual ~LapackQrMemory() {}

    // Matrix
    std::vector<double> mat;

    // The scalar factors of the elementary reflectors
    std::vector<double> tau;

    // qr work array
    std::vector<double> work;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_LAPACK_QR_HPP
