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


#ifndef CASADI_CSPARSE_CHOLESKY_INTERNAL_HPP
#define CASADI_CSPARSE_CHOLESKY_INTERNAL_HPP

/** \defgroup plugin_LinearSolver_csparsecholesky
   * LinearSolver with CSparseCholesky Interface
*/

/** \pluginsection{LinearSolver,csparsecholesky} */

/// \cond INTERNAL
#include <cs.h>
#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/csparse/casadi_linearsolver_csparsecholesky_export.h>

namespace casadi {

  /** \brief \pluginbrief{LinearSolver,csparsecholesky}
   *
   *
   @copydoc LinearSolver_doc
   @copydoc plugin_LinearSolver_csparsecholesky
   *
   */
  class CASADI_LINEARSOLVER_CSPARSECHOLESKY_EXPORT
  CSparseCholeskyInternal : public LinearSolverInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    CSparseCholeskyInternal(const std::string& name, const Sparsity& sp, int nrhs);

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const std::string& name, const Sparsity& sp, int nrhs) {
      return new CSparseCholeskyInternal(name, sp, nrhs);
    }

    // Destructor
    virtual ~CSparseCholeskyInternal();

    // Initialize the solver
    virtual void init();

    // Factorize the matrix
    virtual void linsol_prepare();

    // Solve the system of equations
    virtual void linsol_solve(double* x, int nrhs, bool tr);

    // Solve the system of equations <tt>Lx = b</tt>
    virtual void linsol_solveL(double* x, int nrhs, bool tr);

    /// Obtain a symbolic Cholesky factorization
    virtual Sparsity linsol_cholesky_sparsity(bool tr) const;

    /// Obtain a numeric Cholesky factorization
    virtual DMatrix linsol_cholesky(bool tr) const;

    // The transpose of linear system in form (CCS)
    cs AT_;

    // The symbolic factorization
    css *S_;

    // The numeric factorization
    csn *L_;

    // Temporary
    std::vector<double> temp_;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    virtual const char* plugin_name() const { return "csparsecholesky";}
  };

} // namespace casadi

/// \endcond
#endif // CASADI_CSPARSE_CHOLESKY_INTERNAL_HPP
