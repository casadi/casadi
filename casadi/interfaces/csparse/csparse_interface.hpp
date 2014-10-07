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


#ifndef CASADI_CSPARSE_INTERFACE_HPP
#define CASADI_CSPARSE_INTERFACE_HPP

/** \defgroup plugin_LinearSolver_csparse
 * LinearSolver with CSparse Interface
*/

/** \pluginsection{LinearSolver,csparse} */

/// \cond INTERNAL
#include <cs.h>
#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/csparse/casadi_linearsolver_csparse_export.h>

namespace casadi {

  /** \brief \pluginbrief{LinearSolver,csparse}
   * @copydoc LinearSolver_doc
   * @copydoc plugin_LinearSolver_csparse
   */
  class CASADI_LINEARSOLVER_CSPARSE_EXPORT CsparseInterface : public LinearSolverInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    CsparseInterface(const Sparsity& sp, int nrhs);

    // Copy constructor
    CsparseInterface(const CsparseInterface& linsol);

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new CsparseInterface(sp, nrhs);}

    // Destructor
    virtual ~CsparseInterface();

    // Initialize the solver
    virtual void init();

    // Factorize the matrix
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    // Clone
    virtual CsparseInterface* clone() const;

    // Has the solve function been called once
    bool called_once_;

    // The linear system CSparse form (CCS)
    cs A_;

    // The symbolic factorization
    css *S_;

    // The numeric factorization
    csn *N_;

    // Temporary
    std::vector<double> temp_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond

#endif // CASADI_CSPARSE_INTERFACE_HPP
