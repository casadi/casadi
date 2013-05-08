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

#ifndef SQP_INTERNAL_HPP
#define SQP_INTERNAL_HPP

#include "sqp_method.hpp"
#include "symbolic/fx/nlp_solver_internal.hpp"
#include "symbolic/fx/qp_solver.hpp"
#include <deque>

namespace CasADi{
    
class SQPInternal : public NLPSolverInternal{

public:
  explicit SQPInternal(const FX& nlp);
  virtual ~SQPInternal();
  virtual SQPInternal* clone() const{ return new SQPInternal(*this);}
  
  virtual void init();
  virtual void evaluate(int nfdir, int nadir);
  
  /// QP solver for the subproblems
  QPSolver qp_solver_;

  /// Exact Hessian?
  bool exact_hessian_;

  /// maximum number of sqp iterations
  int maxiter_; 

  /// Memory size of L-BFGS method
  int lbfgs_memory_;
  /// Tolerance of primal infeasibility
  double tol_pr_;
  /// Tolerance of dual infeasibility
  double tol_du_;

  /// Linesearch parameters
  //@{
  double sigma_;
  double c1_;
  double beta_;
  int maxiter_ls_;
  int merit_memsize_;
  //@}

  /// Hessian regularization
  double reg_;
  
  /// Access QPSolver
  const QPSolver getQPSolver() const { return qp_solver_;}
  
  /// Lagrange multipliers of the NLP
  std::vector<double> mu_, mu_x_;
  
  /// Current cost function value
  double fk_;
  
  /// Current and previous linearization point and candidate
  std::vector<double> x_, x_old_, x_cand_;
  
  /// Lagrange gradient in the next iterate
  std::vector<double> gLag_, gLag_old_;
  
  /// Constraint function value
  std::vector<double> gk_, gk_cand_;
  
  /// Gradient of the objective function
  std::vector<double> gf_;

  /// BFGS update function
  enum BFGSMdoe{ BFGS_BK, BFGS_X, BFGS_X_OLD, BFGS_GLAG, BFGS_GLAG_OLD, BFGS_NUM_IN}; 
  FX bfgs_;
  
  /// Initial Hessian approximation (BFGS)
  DMatrix B_init_;
  
  /// Current Hessian approximation
  DMatrix Bk_;
  
  // Current Jacobian
  DMatrix Jk_;

  // Bounds of the QP
  std::vector<double> qp_LBA_, qp_UBA_, qp_LBX_, qp_UBX_;

  // QP solution
  std::vector<double> dx_, qp_DUAL_X_, qp_DUAL_A_;

  /// Regularization
  bool regularize_;

  // Storage for merit function
  std::deque<double> merit_mem_;

  /// Print iteration header
  void printIteration(std::ostream &stream);
  
  /// Print iteration
  void printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, 
                      double dx_norm, double reg, int ls_trials, bool ls_success);

  // Reset the Hessian or Hessian approximation
  void reset_h();

  // Evaluate the gradient of the objective
  virtual void eval_f(const std::vector<double>& x, double& f);
  
  // Evaluate the gradient of the objective
  virtual void eval_grad_f(const std::vector<double>& x, double& f, std::vector<double>& grad_f);
  
  // Evaluate the constraints
  virtual void eval_g(const std::vector<double>& x, std::vector<double>& g);

  // Evaluate the Jacobian of the constraints
  virtual void eval_jac_g(const std::vector<double>& x, std::vector<double>& g, Matrix<double>& J);

  // Evaluate the Hessian of the Lagrangian
  virtual void eval_h(const std::vector<double>& x, const std::vector<double>& lambda, double sigma, Matrix<double>& H);
  
  // Calculate the regularization parameter using Gershgorin theorem
  double getRegularization(const Matrix<double>& H);
  
  // Regularize by adding a multiple of the identity
  void regularize(Matrix<double>& H, double reg);
  
  // Solve the QP subproblem
  virtual void solve_QP(const Matrix<double>& H, const std::vector<double>& g,
                        const std::vector<double>& lbx, const std::vector<double>& ubx,
                        const Matrix<double>& A, const std::vector<double>& lbA, const std::vector<double>& ubA,
                        std::vector<double>& x_opt, std::vector<double>& lambda_x_opt, std::vector<double>& lambda_A_opt);
  
  // Calculate the L1-norm of the primal infeasibility
  double primalInfeasibility(const std::vector<double>& x, const std::vector<double>& lbx, const std::vector<double>& ubx,
                             const std::vector<double>& g, const std::vector<double>& lbg, const std::vector<double>& ubg);
  
  /// Calculates inner_prod(x,mul(A,x))
  static double quad_form(const std::vector<double>& x, const DMatrix& A);
  
};

} // namespace CasADi

#endif //SQP_INTERNAL_HPP
