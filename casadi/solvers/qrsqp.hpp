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


#ifndef CASADI_QRSQP_HPP
#define CASADI_QRSQP_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/solvers/casadi_nlpsol_qrsqp_export.h>

/** \defgroup plugin_Nlpsol_qrsqp
 A textbook SQPMethod
*/

/** \pluginsection{Nlpsol,qrsqp} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_QRSQP_EXPORT QrsqpMemory : public NlpsolMemory {
    /// Current and previous linearization point and candidate
    double *z_cand;

    /// Lagrange gradient in the next iterate
    double *gLag, *gLag_old;

    /// Gradient of the objective function
    double *gf;

    // Bounds of the QP
    double *lbdz, *ubdz;

    // QP solution
    double *dz, *dlam;

    // Current Jacobian
    double *Jk;

    /// Current Hessian approximation
    double *Bk;

    /// Hessian regularization
    double reg;

    /// Linesearch parameters
    double sigma;

    // Storage for merit function
    double* merit_mem;
    size_t merit_ind;

    /// Last return status
    const char* return_status;

    /// Iteration count
    int iter_count;
  };

  /** \brief  \pluginbrief{Nlpsol,sqsqp}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_qrsqp
  */
  class CASADI_NLPSOL_QRSQP_EXPORT Qrsqp : public Nlpsol {
  public:
    explicit Qrsqp(const std::string& name, const Function& nlp);
    ~Qrsqp() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "qrsqp";}

    // Name of the class
    std::string class_name() const override { return "Qrsqp";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Qrsqp(name, nlp);
    }

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new QrsqpMemory();}

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<QrsqpMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    /// QP solver for the subproblems
    Function qpsol_;

    /// Exact Hessian?
    bool exact_hessian_;

    /// Maximum, minimum number of SQP iterations
    casadi_int max_iter_, min_iter_;

    /// Memory size of L-BFGS method
    casadi_int lbfgs_memory_;

    /// Tolerance of primal and dual infeasibility
    double tol_pr_, tol_du_;

    /// Minimum step size allowed
    double min_step_size_;

    /// Linesearch parameters
    ///@{
    double c1_;
    double beta_;
    casadi_int max_iter_ls_;
    casadi_int merit_memsize_;
    ///@}

    // Print options
    bool print_header_, print_iteration_;

    // Hessian sparsity
    Sparsity Hsp_;

    // Jacobian sparsity
    Sparsity Asp_;

    /// Regularization
    bool regularize_;

    /// Access Conic
    const Function getConic() const { return qpsol_;}

    /// Print iteration header
    void print_iteration() const;

    /// Print iteration
    void print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf,
                         double dx_norm, double rg, casadi_int ls_trials, bool ls_success) const;

    // Solve the QP subproblem
    virtual void solve_QP(QrsqpMemory* m, const double* H, const double* g,
                          const double* lbdz, const double* ubdz,
                          const double* A, double* x_opt, double* dlam) const;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_QRSQP_HPP
