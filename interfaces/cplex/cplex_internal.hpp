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
#ifndef CPLEX_INTERNAL_HPP
#define CPLEX_INTERNAL_HPP

#include "ilcplex/cplex.h"
#include "symbolic/function/qp_solver_internal.hpp"

#include <string>

/// \cond INTERNAL

namespace CasADi{

  /** Internal class for CplexSolver
      @copydoc QPSolver_doc
  */
  class CplexInternal : public QPSolverInternal{
    friend class CplexSolver;
  public:
    /** \brief Default constructor */
    explicit CplexInternal();
  
    /// Clone
    virtual CplexInternal* clone() const;
  
    /// Constructor using sparsity patterns
    explicit CplexInternal(const std::vector<Sparsity>& st);

    /// Destructor
    virtual ~CplexInternal();
  
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
     * 6 -> Concurent
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
  
    /// Determines relation >,<,= in the lin. constraints
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

    /// CPLEX-environment
    CPXENVptr env_;
    CPXLPptr lp_;

  };
} // end namespace CasADi
/// \endcond
#endif //CPLEX_INTERNAL_HPP
