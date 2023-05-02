/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *    Copyright (C) 2022-2023 David Kiessling
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


#ifndef CASADI_FEASIBLESQPMETHOD_HPP
#define CASADI_FEASIBLESQPMETHOD_HPP

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/solvers/casadi_nlpsol_feasiblesqpmethod_export.h>

/** \defgroup plugin_Nlpsol_feasiblesqpmethod
    \par
 A textbook FeasibleSQPMethod

    \identifier{241} */

/** \pluginsection{Nlpsol,feasiblesqpmethod} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_FEASIBLESQPMETHOD_EXPORT FeasiblesqpmethodMemory : public NlpsolMemory {
    // Problem data structure
    casadi_feasiblesqpmethod_data<double> d;
    /// Hessian regularization
    double reg;

    /// Linesearch parameters
    double sigma;

    casadi_int merit_ind;
    /// Last return status
    const char* return_status;

    /// Iteration count
    int iter_count;
  };

  /** \brief  \pluginbrief{Nlpsol,feasiblesqpmethod}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_feasiblesqpmethod
  */
  class CASADI_NLPSOL_FEASIBLESQPMETHOD_EXPORT Feasiblesqpmethod : public Nlpsol {
  public:
    explicit Feasiblesqpmethod(const std::string& name, const Function& nlp);
    ~Feasiblesqpmethod() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "feasiblesqpmethod";}

    // Name of the class
    std::string class_name() const override { return "Feasiblesqpmethod";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Feasiblesqpmethod(name, nlp);
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
    void* alloc_mem() const override { return new FeasiblesqpmethodMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<FeasiblesqpmethodMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    double eval_m_k(void* mem) const;

    double eval_tr_ratio(double val_f, double val_f_corr, double val_m_k) const;

    void tr_update(void* mem, double& tr_rad, double tr_ratio) const;

    int step_update(void* mem, double tr_ratio) const;

    // function to get feasible iterate
    int feasibility_iterations(void* mem, double tr_rad) const;

    void anderson_acc_step_update(void* mem, casadi_int iter_index) const;

    void anderson_acc_init_memory(void* mem, double* step, double* iterate) const;

    void anderson_acc_update_memory(void* mem, double* step, double* iterate) const;

    // Solve the NLP
    int solve(void* mem) const override;

    // Memory structure
    casadi_feasiblesqpmethod_prob<double> p_;

    /// QP solver for the subproblems
    Function qpsol_;

    /// QP solver for elastic mode subproblems
    Function qpsol_ela_;

    /// Exact Hessian?
    bool exact_hessian_;

    /// Exact Hessian?
    bool use_sqp_;

    /// Use Anderson Acceleration
    bool use_anderson_;

    /// Maximum block size of Hessian
    casadi_int block_size_ = 0;

    /// Maximum, minimum number of SQP iterations
    casadi_int max_iter_, min_iter_;

    /// Memory size of L-BFGS method
    casadi_int lbfgs_memory_;

    // Memory size of Anderson acceleration
    casadi_int sz_anderson_memory_;

    /// Tolerance of primal and dual infeasibility
    double tol_pr_, tol_du_;

    /// Initialize feasible qp's
    bool init_feasible_;

    /// Linesearch parameters
    ///@{
    // double c1_;
    // double beta_;
    // casadi_int max_iter_ls_;
    // casadi_int merit_memsize_;
    ///@}

    // Print options
    bool print_header_, print_iteration_, print_status_;

    // Hessian Sparsity
    Sparsity Hsp_;

    // Jacobian sparsity
    Sparsity Asp_;

    /// Data for convexification
    ConvexifyData convexify_data_;

    /// convexify?
    bool convexify_;

    // -------- FROM HERE OPTIONS FOR FP-SQP ------------------

    // tolerances
    double optim_tol_, feas_tol_;

    // trust-region parameters
    double tr_eta1_, tr_eta2_;
    double tr_alpha1_, tr_alpha2_;
    double tr_tol_;
    double tr_rad0_;
    double tr_acceptance_;
    double tr_rad_min_, tr_rad_max_;

    std::vector<double> tr_scale_vector_;

    // inner iterations
    double contraction_acceptance_value_;
    casadi_int watchdog_;
    casadi_int max_inner_iter_;


    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /// Access Conic
    const Function getConic() const { return qpsol_;}

    /// Print iteration header
    void print_iteration() const;

    /// Print iteration
    void print_iteration(casadi_int iter, double obj, double m_k,
                         double tr_ratio, double pr_inf, double du_inf,
                         double dx_norm, double rg, double tr_rad,
                         std::string info) const;

    // Solve the QP subproblem: mode 0 = normal, mode 1 = SOC
    virtual int solve_QP(FeasiblesqpmethodMemory* m, const double* H, const double* g,
                          const double* lbdz, const double* ubdz,
                          const double* A,
                          double* x_opt, double* dlam, int mode) const;

    // Solve the LP
    virtual int solve_LP(FeasiblesqpmethodMemory* m, const double* g,
                          const double* lbdz, const double* ubdz,
                          const double* A,
                          double* x_opt, double* dlam, int mode) const;

    // Solve the QP subproblem
    void codegen_qp_solve(CodeGenerator& cg, const std::string& H, const std::string& g,
              const std::string& lbdz, const std::string& ubdz,
              const std::string& A, const std::string& x_opt,
              const std::string& dlam, int mode) const;

    void codegen_tr_update(CodeGenerator& cg,
      const std::string& tr_rad, const std::string& tr_ratio) const;

    void codegen_eval_m_k(CodeGenerator& cg) const;

    void codegen_eval_tr_ratio(CodeGenerator& cg,
      const std::string& val_f, const std::string& val_f_corr, const std::string& val_m_k) const;

    void codegen_step_update(CodeGenerator& cg, const std::string& tr_ratio) const;

    void codegen_feasibility_iterations(CodeGenerator& cg, const std::string& tr_rad) const;

    // Solve the QP subproblem
    // void codegen_qp_ela_solve(CodeGenerator& cg, const std::string& H, const std::string& g,
    //           const std::string& lbdz, const std::string& ubdz,
    //           const std::string& A, const std::string& x_opt, const std::string& dlam) const;

    // Execute elastic mode: mode 0 = normal, mode 1 = SOC
    // void codegen_solve_elastic_mode(CodeGenerator& cg, int mode) const;

    // Codegen to calculate gama_1
    // void codegen_calc_gamma_1(CodeGenerator& cg) const;


    // Calculate gamma_1
    // double calc_gamma_1(FeasiblesqpmethodMemory* m) const;

    /// A documentation string
    static const std::string meta_doc;


    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Feasiblesqpmethod(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit Feasiblesqpmethod(DeserializingStream& s);

  private:
    void set_feasiblesqpmethod_prob();
  };

} // namespace casadi
/// \endcond
#endif // CASADI_FEASIBLESQPMETHOD_HPP
