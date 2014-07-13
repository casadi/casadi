/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef CASADI_CSPARSE_INTERNAL_HPP
#define CASADI_CSPARSE_INTERNAL_HPP

/// \cond INTERNAL
extern "C" {
#include "external_packages/CSparse/Include/cs.h"
}
#include "casadi/core/function/linear_solver_internal.hpp"
#include <casadi/interfaces/csparse/casadi_linearsolver_csparse_export.h>

namespace casadi {

  /** \brief  LinearSolver with CSparse Interface
   *
   @copydoc LinearSolver_doc
   *
   * CSparse is an casadi::Function mapping from 2 inputs
   * [ A (matrix), b (vector)] to one output [x (vector)].
   *
   * The usual procedure to use CSparse is: \n
   *  -# init()
   *  -# set the first input (A)
   *  -# prepare()
   *  -# set the second input (b)
   *  -# solve()
   *  -# Repeat steps 4 and 5 to work with other b vectors.
   *
   * The method evaluate() combines the prepare() and solve()
   * step and is therefore more expensive if A is invariant.
   *
   */
  class CASADI_LINEARSOLVER_CSPARSE_EXPORT CSparseInternal : public LinearSolverInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    CSparseInternal(const Sparsity& sp, int nrhs);

    // Copy constructor
    CSparseInternal(const CSparseInternal& linsol);

    /** \brief  Create a new LinearSolver */
    static LinearSolverInternal* creator(const Sparsity& sp, int nrhs)
    { return new CSparseInternal(sp, nrhs);}

    // Destructor
    virtual ~CSparseInternal();

    // Initialize the solver
    virtual void init();

    // Factorize the matrix
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    // Clone
    virtual CSparseInternal* clone() const;

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


  };

} // namespace casadi

/// \endcond

#endif // CASADI_CSPARSE_INTERNAL_HPP
