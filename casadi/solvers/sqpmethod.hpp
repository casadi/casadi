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

#include "casadi/core/function/nlpsol_impl.hpp"
#include <deque>

#include <casadi/solvers/casadi_nlpsol_sqpmethod_export.h>

/** \defgroup plugin_Nlpsol_sqpmethod
 A textbook SQPMethod
*/

/** \pluginsection{Nlpsol,sqpmethod} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_SQPMETHOD_EXPORT SqpmethodMemory : public NlpsolMemory {
    /// Current cost function value
    double fk;

    /// Lagrange multipliers of the NLP
    double *mu, *mu_x;

    /// Current and previous linearization point and candidate
    double *xk, *x_old, *x_cand;

    /// Lagrange gradient in the next iterate
    double *gLag, *gLag_old;

    /// Constraint function value
    double *gk, *gk_cand;

    /// Gradient of the objective function
    double *gf;

    // Bounds of the QP
    double *qp_LBA, *qp_UBA, *qp_LBX, *qp_UBX;

    // QP solution
    double *dx, *qp_DUAL_X, *qp_DUAL_A;

    // Current Jacobian
    double *Jk;

    /// Current Hessian approximation
    double *Bk;

    /// Hessian regularization
    double reg;

    /// Linesearch parameters
    double sigma;

    // Storage for merit function
    std::deque<double> merit_mem;

    /// Last return status
    const char* return_status;

    // Accumulated time since last reset:
    double t_eval_f; // time spent in eval_f
    double t_eval_grad_f; // time spent in eval_grad_f
    double t_eval_g; // time spent in eval_g
    double t_eval_jac_g; // time spent in eval_jac_g
    double t_eval_h; // time spent in eval_h
    double t_callback_fun;  // time spent in callback function
    double t_callback_prepare; // time spent in callback preparation
    double t_mainloop; // time spent in the main loop of the solver

    // Accumulated counts since last reset:
    int n_eval_f; // number of calls to eval_f
    int n_eval_grad_f; // number of calls to eval_grad_f
    int n_eval_g; // number of calls to eval_g
    int n_eval_jac_g; // number of calls to eval_jac_g
    int n_eval_h; // number of calls to eval_h
  };

  /** \brief  \pluginbrief{Nlpsol,sqpmethod}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_sqpmethod
  */
  class CASADI_NLPSOL_SQPMETHOD_EXPORT Sqpmethod : public Nlpsol {
  public:
    explicit Sqpmethod(const std::string& name, const XProblem& nlp);
    virtual ~Sqpmethod();

  // Get name of the plugin
  virtual const char* plugin_name() const { return "sqpmethod";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const XProblem& nlp) {
      return new Sqpmethod(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual Memory* memory() const { return new SqpmethodMemory();}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(Memory& mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Solve the NLP
    virtual void solve(Memory& mem) const;

    /// QP solver for the subproblems
    Function qpsol_;

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
    double c1_;
    double beta_;
    int max_iter_ls_;
    int merit_memsize_;
    ///@}

    // Print options
    bool print_header_;
    bool print_time_;

    /// BFGS update function
    enum BFGSMdoe { BFGS_BK, BFGS_X, BFGS_X_OLD, BFGS_GLAG, BFGS_GLAG_OLD, BFGS_NUM_IN};
    Function bfgs_;

    // Hessian sparsity
    Sparsity Hsp_;

    // Jacobian sparsity
    Sparsity Asp_;

    /// Initial Hessian approximation (BFGS)
    DM B_init_;

    /// Regularization
    bool regularize_;

    /// Access Qpsol
    const Function getQpsol() const { return qpsol_;}

    /// Print iteration header
    void printIteration(std::ostream &stream) const;

    /// Print iteration
    void printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf,
                        double dx_norm, double reg, int ls_trials, bool ls_success) const;

    // Reset the Hessian or Hessian approximation
    void reset_h(SqpmethodMemory& m) const;

    // Evaluate the gradient of the objective
    virtual double eval_f(SqpmethodMemory& m, const double* x) const;

    // Evaluate the gradient of the objective
    virtual void eval_grad_f(SqpmethodMemory& m, const double* x, double* f,
                             double* grad_f) const;

    // Evaluate the constraints
    virtual void eval_g(SqpmethodMemory& m, const double* x, double* g) const;

    // Evaluate the Jacobian of the constraints
    virtual void eval_jac_g(SqpmethodMemory& m, const double* x, double* g, double* J) const;

    // Evaluate the Hessian of the Lagrangian
    virtual void eval_h(SqpmethodMemory& m, const double* x, const double* lambda,
                        double sigma, double* H) const;

    // Calculate the regularization parameter using Gershgorin theorem
    double getRegularization(const double* H) const;

    // Regularize by adding a multiple of the identity
    void regularize(double* H, double reg) const;

    // Solve the QP subproblem
    virtual void solve_QP(SqpmethodMemory& m, const double* H, const double* g,
                          const double* lbx, const double* ubx,
                          const double* A, const double* lbA, const double* ubA,
                          double* x_opt, double* lambda_x_opt, double* lambda_A_opt) const;

    // Calculate the L1-norm of the primal infeasibility
    double primalInfeasibility(const double* x,
                               const double* lbx, const double* ubx,
                               const double* g, const double* lbg, const double* ubg) const;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SQPMETHOD_HPP
