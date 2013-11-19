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
#include "symbolic/fx/qp_solver_internal.hpp"

#include <string>

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
  explicit CplexInternal(const std::vector<CRSSparsity>& st);

  /// Destructor
  virtual ~CplexInternal();
  
  /// Free Cplex memory
  void freeCplex();

  // Initialize the solver
  virtual void init();

  // Solve the QP
  virtual void evaluate();

  /// Converts CasADi sparsity to Cplex sparsity
  static void toCplexSparsity(const CRSSparsity& sp_trans, std::vector<int> &matbeg, std::vector<int>& matcnt, std::vector<int>& matind);
  
protected:
  
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
  /// Number of variables
  int NUMCOLS_;
  /// Number of constrains (altogether)
  int NUMROWS_;
  /// Nature of problem (always minimization)
  int objsen_;
  
  /// Linear term of objective
  double* obj_;

  /// Determines relation >,<,= in the lin. constraints
  std::vector<char> sense_;

  /// Coefficients of matrix A (constraint Jacobian)
  std::vector<int> matbeg_;
  std::vector<int> matcnt_;
  std::vector<int> matind_;
  std::vector<double> matval_;

  /// Transposed sparsity pattern mapping from A to AT
  std::vector<int> AT_nonzero_mapping_;

  /// Right-hand side of constraints
  std::vector<double> rhs_;
  /// Range of constraints
  std::vector<double> rngval_;

  /// Simple bounds on x
  double  *lb_;
  double  *ub_;
  
  /// Coefficients of matrix H (objective Hessian)
  std::vector<int> qmatbeg_;
  std::vector<int> qmatcnt_;
  std::vector<int> qmatind_;
  
  /// Storage for basis info of primal variables
  std::vector<int> cstat_;
  /// Storage for basis info of slack variables
  std::vector<int> rstat_;

  /// CPLEX-environment
  CPXENVptr env_;
  CPXLPptr lp_;

};
} // end namespace CasADi

#endif //CPLEX_INTERNAL_HPP
