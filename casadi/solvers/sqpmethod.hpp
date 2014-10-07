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


#ifndef CASADI_SQPMETHOD_HPP
#define CASADI_SQPMETHOD_HPP

#include "casadi/core/function/nlp_solver_internal.hpp"
#include "casadi/core/function/qp_solver.hpp"
#include <deque>

#include <casadi/solvers/casadi_nlpsolver_sqpmethod_export.h>

/** \defgroup plugin_NlpSolver_sqpmethod
 A textbook SQPMethod
*/

/** \pluginsection{NlpSolver,sqpmethod} */

/// \cond INTERNAL
namespace casadi {

  /** \brief  \pluginbrief{NlpSolver,sqpmethod}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_NlpSolver_sqpmethod
  */
  class CASADI_NLPSOLVER_SQPMETHOD_EXPORT Sqpmethod : public NlpSolverInternal {
  public:
    explicit Sqpmethod(const Function& nlp);
    virtual ~Sqpmethod();
    virtual Sqpmethod* clone() const { return new Sqpmethod(*this);}

    /** \brief  Create a new NLP Solver */
    static NlpSolverInternal* creator(const Function& nlp)
    { return new Sqpmethod(nlp);}

    virtual void init();
    virtual void evaluate();

    /// QP solver for the subproblems
    QpSolver qp_solver_;

    /// Exact Hessian?
    bool exact_hessian_;

    /// maximum number of sqp iterations
    int max_iter_;

    /// Memory size of L-BFGS method
    int lbfgs_memory_;
    /// Tolerance of primal infeasibility
    double tol_pr_;
    /// Tolerance of dual infeasibility
    double tol_du_;

    /// Minimum step size allowed
    double min_step_size_;

    /// Linesearch parameters
    ///@{
    double sigma_;
    double c1_;
    double beta_;
    int max_iter_ls_;
    int merit_memsize_;
    ///@}

    /// Hessian regularization
    double reg_;

    /// Access QpSolver
    const QpSolver getQpSolver() const { return qp_solver_;}

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
    enum BFGSMdoe { BFGS_BK, BFGS_X, BFGS_X_OLD, BFGS_GLAG, BFGS_GLAG_OLD, BFGS_NUM_IN};
    Function bfgs_;

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
    virtual void eval_jac_g(const std::vector<double>& x, std::vector<double>& g,
      Matrix<double>& J);

    // Evaluate the Hessian of the Lagrangian
    virtual void eval_h(const std::vector<double>& x, const std::vector<double>& lambda,
                        double sigma, Matrix<double>& H);

    // Calculate the regularization parameter using Gershgorin theorem
    double getRegularization(const Matrix<double>& H);

    // Regularize by adding a multiple of the identity
    void regularize(Matrix<double>& H, double reg);

    // Solve the QP subproblem
    virtual void solve_QP(const Matrix<double>& H, const std::vector<double>& g,
                          const std::vector<double>& lbx, const std::vector<double>& ubx,
                          const Matrix<double>& A, const std::vector<double>& lbA,
                          const std::vector<double>& ubA,
                          std::vector<double>& x_opt, std::vector<double>& lambda_x_opt,
                          std::vector<double>& lambda_A_opt);

    // Calculate the L1-norm of the primal infeasibility
    double primalInfeasibility(const std::vector<double>& x, const std::vector<double>& lbx,
                               const std::vector<double>& ubx,
                               const std::vector<double>& g, const std::vector<double>& lbg,
                               const std::vector<double>& ubg);

    /// Calculates <tt>inner_prod(x, mul(A, x))</tt>
    static double quad_form(const std::vector<double>& x, const DMatrix& A);

    // Accumulated time since last reset:
    double t_eval_f_; // time spent in eval_f
    double t_eval_grad_f_; // time spent in eval_grad_f
    double t_eval_g_; // time spent in eval_g
    double t_eval_jac_g_; // time spent in eval_jac_g
    double t_eval_h_; // time spent in eval_h
    double t_callback_fun_;  // time spent in callback function
    double t_callback_prepare_; // time spent in callback preparation
    double t_mainloop_; // time spent in the main loop of the solver

    // Accumulated counts since last reset:
    int n_eval_f_; // number of calls to eval_f
    int n_eval_grad_f_; // number of calls to eval_grad_f
    int n_eval_g_; // number of calls to eval_g
    int n_eval_jac_g_; // number of calls to eval_jac_g
    int n_eval_h_; // number of calls to eval_h

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SQPMETHOD_HPP
