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
  explicit SCPgenInternal(const FX& F, const FX& G, const FX& H, const FX& J);
  virtual ~SCPgenInternal();
  virtual SCPgenInternal* clone() const{ return new SCPgenInternal(*this);}
  
  virtual void init();
  virtual void evaluate(int nfdir, int nadir);

  // Codegen function
  void dynamicCompilation(FX& f, FX& f_gen, std::string fname, std::string fdescr);

  // Calculate the L1-norm of the primal infeasibility
  double primalInfeasibility();

  // Calculate the L1-norm of the dual infeasibility
  double dualInfeasibility();

  // Print iteration header
  void printIteration(std::ostream &stream);
  
  // Print iteration
  void printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, 
                      double reg, int ls_trials, bool ls_success);

  // Evaluate the residual function
  void eval_res();

  // Form the condensed QP
  void eval_qpf();

  // Regularize the condensed QP
  void regularize();

  // Solve the QP to get the (full) step
  void solve_qp();

  // Perform the line-search to take the step
  void line_search(int& ls_iter, bool& ls_success);

  // Evaluate the step expansion
  void eval_exp();



  
  /// QP solver for the subproblems
  QPSolver qp_solver_;

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
  double c1_;
  double beta_;
  int maxiter_ls_;
  int merit_memsize_;
  //@}

  /// Enable Code generation
  bool codegen_;

  /// Access QPSolver
  const QPSolver getQPSolver() const { return qp_solver_;}
  
  /// Regularization
  bool regularize_;

  // Storage for merit function
  std::deque<double> merit_mem_;

  // Options
  double reg_threshold_;

  /// stopping criterion for the stepsize
  double toldx_;
  
  /// stopping criterion for the lagrangian gradient
  double tolgl_;
  
  /// Generate initial guess for lifted variables
  FX vinit_fcn_;

  /// Residual function
  FX res_fcn_;
 
  /// Quadratic approximation
  FX qp_fcn_;

  /// Step expansion
  FX exp_fcn_;
  
  /// Dimensions
  int ngL_;

  // Objective value
  double obj_k_;

  // Simple and nonlinear bounds
  std::vector<double> lbu_, ubu_, g_, lbg_, ubg_, gL_;

  /// Multipliers for the nonlinear bounds
  std::vector<double> lambda_g_, dlambda_g_;

  int res_lam_g_;
  int res_obj_, res_gl_, res_g_;

  int z_p_, z_lam_g_;
  int z_obj_, z_gl_, z_g_;

  int qpf_p_, qpf_lam_g_;
  int qpf_b_obj_, qpf_B_obj_, qpf_b_g_, qpf_B_g_;

  int exp_p_, exp_du_, exp_dlam_g_, exp_lam_g_;
  int exp_osens_, exp_curve_;

  struct Var{
    int n;
    MX v;

    int res_var, res_lam;
    int res_d, res_lam_d;

    int z_var, z_lam;
    int z_def, z_defL;

    int qpf_var, qpf_lam;
    int qpf_def, qpf_defL;

    int exp_var, exp_lam;
    int exp_def, exp_defL;

    std::vector<double> step, init, opt, lam, dlam;
    std::vector<double> res, resL;
    
  };
  std::vector<Var> x_;  

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

  // QP
  DMatrix qpH_, qpA_;
  std::vector<double> qpG_, qpB_;
};

} // namespace CasADi

#endif //SCPGEN_INTERNAL_HPP
