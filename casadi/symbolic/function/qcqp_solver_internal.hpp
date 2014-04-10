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

#ifndef QCQP_SOLVER_INTERNAL_HPP
#define QCQP_SOLVER_INTERNAL_HPP

#include "qcqp_solver.hpp"
#include "function_internal.hpp"


/// \cond INTERNAL
namespace casadi{

/// Internal class
class CASADI_SYMBOLIC_EXPORT QCQPSolverInternal : public FunctionInternal{
  public:

    // Constructor
    QCQPSolverInternal(const std::vector<Sparsity> &st);

    // Destructor
    virtual ~QCQPSolverInternal() = 0;

    // Initialize
    virtual void init();

    // Solve the system of equations
    virtual void evaluate();

    // Solve the system of equations
    virtual void solve();

    /// Set options that make the QP solver more suitable for solving LPs
    virtual void setQPOptions() {}

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void checkInputs() const;

  protected:

    /// Problem structure
    std::vector<Sparsity> st_;

    /// Number of decision variables
    int n_;

    /// The number of linear constraints (counting both equality and inequality) == A.size1()
    int nc_;

    /// The number of quadratic constraints
    int nq_;
};


} // namespace casadi

/// \endcond
#endif //QCQP_SOLVER_INTERNAL_HPP

