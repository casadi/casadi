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

#ifndef CASADI_CPLEX_INTERFACE_HPP
#define CASADI_CPLEX_INTERFACE_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include <casadi/interfaces/cplex/casadi_qpsolver_cplex_export.h>
#include "ilcplex/cplex.h"

#include <string>

/** \defgroup plugin_QpSolver_cplex

      Interface to Cplex solver for sparse Quadratic Programs
*/

/** \pluginsection{QpSolver,cplex} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{QpSolver,cplex}

      @copydoc QpSolver_doc
      @copydoc plugin_QpSolver_cplex

      \author Attila Kozma, Joel Andersson
      \date 2012
  */
  class CASADI_QPSOLVER_CPLEX_EXPORT CplexInterface : public QpSolverInternal {
  public:
    /** \brief Default constructor */
    explicit CplexInterface();

    /// Clone
    virtual CplexInterface* clone() const;

    /** \brief  Create a new QP Solver */
    static QpSolverInternal* creator(const QPStructure& st)
    { return new CplexInterface(st);}

    /// Constructor using sparsity patterns
    explicit CplexInterface(const std::vector<Sparsity>& st);

    /// Destructor
    virtual ~CplexInterface();

    /// Free Cplex memory
    void freeCplex();

    // Initialize the solver
    virtual void init();

    // Solve the QP
    virtual void evaluate();

    // OPTIONS
    /** Which algorithm to use
     * 0 -> Automatic (default)
     * 1 -> Primal simplex
     * 2 -> Dual simplex
     * 3 -> Network optimizer
     * 4 -> Barrier
     * 5 -> Sifting
     * 6 -> Concurrent
     * 7 -> Crossover
     */
    /// Stores which QP algorithm to use
    int qp_method_;

    /// Print to file (for later use)
    bool dump_to_file_;

    /// Indicates if we have to warm-start
    bool is_warm_;

    /// Accuracy
    double tol_;

    /// Nature of problem (always minimization)
    int objsen_;

    /// Determines relation >,<, = in the linear constraints
    std::vector<char> sense_;

    /// Coefficients of matrix A (constraint Jacobian)
    std::vector<int> matcnt_;

    /// Right-hand side of constraints
    std::vector<double> rhs_;

    /// Range of constraints
    std::vector<double> rngval_;

    /// Coefficients of matrix H (objective Hessian)
    std::vector<int> qmatcnt_;

    /// Storage for basis info of primal variables
    std::vector<int> cstat_;

    /// Storage for basis info of slack variables
    std::vector<int> rstat_;

    /// CPLEX environment
    CPXENVptr env_;
    CPXLPptr lp_;

    /// A documentation string
    static const std::string meta_doc;

  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CPLEX_INTERFACE_HPP
