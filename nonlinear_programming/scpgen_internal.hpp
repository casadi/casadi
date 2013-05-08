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

#ifndef SCPGEN_INTERNAL_HPP
#define SCPGEN_INTERNAL_HPP

#include "scpgen.hpp"
#include "symbolic/fx/nlp_solver_internal.hpp"
#include "symbolic/fx/qp_solver.hpp"
#include <deque>

namespace CasADi{
    
class SCPgenInternal : public NLPSolverInternal{

public:
  explicit SCPgenInternal(const FX& nlp);
  virtual ~SCPgenInternal();
  virtual SCPgenInternal* clone() const{ return new SCPgenInternal(*this);}
  
  virtual void init();
  virtual void evaluate(int nfdir, int nadir);

  // Calculate the L1-norm of the primal infeasibility
  double primalInfeasibility();

  // Calculate the L1-norm of the dual infeasibility
  double dualInfeasibility();

  // Print iteration header
  void printIteration(std::ostream &stream);
  
  // Print iteration
  void printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, 
                      double reg, int ls_trials, bool ls_success);

  // Evaluate the matrices in the condensed QP
  void eval_mat();

  // Evaluate the vectors in the condensed QP
  void eval_vec();

  // Evaluate the residual function
  void eval_res();

  // Regularize the condensed QP
  void regularize();

  // Solve the QP to get the (full) step
  void solve_qp();

  // Perform the line-search to take the step
  void line_search(int& ls_iter, bool& ls_success);

  // Evaluate the step expansion
  void eval_exp();
  
  // Timings
  double t_eval_mat_, t_eval_res_, t_eval_vec_, t_eval_exp_, t_solve_qp_, t_mainloop_;
  
  /// QP solver for the subproblems
  QPSolver qp_solver_;

  /// use Gauss-Newton Hessian
  bool gauss_newton_; 

  /// maximum number of sqp iterations
  int maxiter_; 

  /// Memory size of L-BFGS method
  int lbfgs_memory_;

  /// Tolerance on primal infeasibility
  double tol_pr_;

  /// Tolerance on dual infeasibility
  double tol_du_;

  /// Tolerance on regularization
  double tol_reg_;

  /// stopping criterion for the stepsize
  double tol_pr_step_;
  
  /// stopping criterion for the lagrangian gradient
  double tol_gl_;

  /// Linesearch parameters
  //@{
  double c1_;
  double beta_;
  int maxiter_ls_;
  std::vector<double> merit_mem_;
  int merit_memsize_;
  double merit_start_;
  int merit_ind_;
  //@}

  /// Enable Code generation
  bool codegen_;

  /// Access QPSolver
  const QPSolver getQPSolver() const { return qp_solver_;}
  
  /// Regularization
  bool regularize_;

  // Options
  double reg_threshold_;
  
  /// Print timers
  bool print_time_;

  /// Generate initial guess for lifted variables
  FX vinit_fcn_;

  /// Residual function
  FX res_fcn_;
 
  // Function to calculate the matrices in the reduced QP
  FX mat_fcn_;
  int mat_jac_, mat_hes_;

  /// Quadratic approximation
  FX vec_fcn_;
  int vec_gf_, vec_g_;

  /// Step expansion
  FX exp_fcn_;
  
  // Objective value
  double f_;
  std::vector<double> gf_, gL_, b_gn_;

  // Nonlifted variables with bound
  std::vector<double> x_lb_, x_ub_, x_init_, x_opt_, x_step_, x_lam_, x_dlam_;
  MX x_;

  // Parameter
  MX p_;

  // Nonlinear bounds
  std::vector<double> g_, g_lb_, g_ub_, g_lam_, g_dlam_;

  // Residual function io indices
  int res_x_, res_p_, res_g_lam_, res_p_lam_, res_p_d_;
  int res_f_, res_gl_, res_g_;

  // Modifier function io indices
  int mod_x_, mod_p_, mod_g_lam_;
  int mod_f_, mod_gl_, mod_g_;
  int mod_du_, mod_dlam_g_;

  struct Var{
    int n;
    MX v, v_def, v_lam, v_defL;
    MX d, d_def, d_lam, d_defL;

    // Indices of function inputs and outputs
    int res_var, res_lam, res_d, res_lam_d;
    int mod_var, mod_lam, mod_def, mod_defL;
    int exp_def, exp_defL;

    std::vector<double> step, init, opt, lam, dlam;
    std::vector<double> res, resL;
    
  };

  std::vector<Var> v_;

  // Penalty parameter of merit function
  double sigma_;

  // 1-norm of last primal step
  double pr_step_;
  
  // 1-norm of last dual step
  double du_step_;

  // Regularization
  double reg_;
  
  // Names of the components
  std::vector<std::string> name_x_;

  // Components to print
  std::vector<int> print_x_;

  // Message applying to a particular iteration
  std::string iteration_note_;
  
  // QP
  DMatrix qpH_, qpA_;
  std::vector<double> qpB_;

  // Hessian times a step
  std::vector<double> qpH_times_du_;
};

} // namespace CasADi

#endif //SCPGEN_INTERNAL_HPP
