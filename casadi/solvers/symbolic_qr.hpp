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


#ifndef CASADI_SYMBOLIC_QR_HPP
#define CASADI_SYMBOLIC_QR_HPP

#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/solvers/casadi_linearsolver_symbolicqr_export.h>

/** \defgroup plugin_LinearSolver_symbolicqr

       LinearSolver based on QR factorization with sparsity pattern based reordering
      _without_ partial pivoting
*/

/** \pluginsection{LinearSolver,symbolicqr} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{LinearSolver,symbolicqr}

      @copydoc LinearSolver_doc
      @copydoc plugin_LinearSolver_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINEARSOLVER_SYMBOLICQR_EXPORT SymbolicQr
    : public LinearSolverInternal {
  public:
    // Constructor
    SymbolicQr(const Sparsity& sparsity, int nrhs);

    // Destructor
    virtual ~SymbolicQr();

    /** \brief  Clone */
    virtual SymbolicQr* clone() const { return new SymbolicQr(*this);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new SymbolicQr(sp, nrhs);}

    // Initialize
    virtual void init();

    // Prepare the factorization
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(std::ostream &stream, const std::string& type,
                                      CodeGenerator& gen) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(std::ostream &stream, const std::string& type,
                              CodeGenerator& gen) const;

    /** \brief Evaluate symbolically (SX) */
    virtual void evaluateSXGen(const SXPtrV& input, SXPtrV& output, bool tr);

    // Factorization function
    Function fact_fcn_;

    // Solve function
    Function solv_fcn_N_, solv_fcn_T_;

    // Storage for QR factorization
    DMatrix Q_, R_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_SYMBOLIC_QR_HPP

