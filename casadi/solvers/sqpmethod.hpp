/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, Kobe Bergmans
 *                            KU Leuven. All rights reserved.
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

#include "casadi/core/nlpsol_impl.hpp"
#include <casadi/solvers/casadi_nlpsol_sqpmethod_export.h>

/** \defgroup plugin_Nlpsol_sqpmethod Title
    \par

 A textbook SQPMethod

    \identifier{22x} */

/** \pluginsection{Nlpsol,sqpmethod} */

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_SQPMETHOD_EXPORT SqpmethodMemory : public NlpsolMemory {
    // Problem data structure
    casadi_sqpmethod_data<double> d;
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

  /** \brief  \pluginbrief{Nlpsol,sqpmethod}
  *  @copydoc NLPSolver_doc
  *  @copydoc plugin_Nlpsol_sqpmethod
  */
  class CASADI_NLPSOL_SQPMETHOD_EXPORT Sqpmethod : public Nlpsol {
  public:
    explicit Sqpmethod(const std::string& name, const Function& nlp);
    ~Sqpmethod() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "sqpmethod";}

    // Name of the class
    std::string class_name() const override { return "Sqpmethod";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new Sqpmethod(name, nlp);
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
    void* alloc_mem() const override { return new SqpmethodMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<SqpmethodMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    // Solve the NLP
    int solve(void* mem) const override;

    // Memory structure
    casadi_sqpmethod_prob<double> p_;

    /// QP solver for the subproblems
    Function qpsol_;

    /// QP solver for elastic mode subproblems
    Function qpsol_ela_;

    /// Exact Hessian?
    bool exact_hessian_;

    /// Maximum block size of Hessian
    casadi_int block_size_ = 0;

    /// Maximum, minimum number of SQP iterations
    casadi_int max_iter_, min_iter_;

    /// Memory size of L-BFGS method
    casadi_int lbfgs_memory_;

    /// Tolerance of primal and dual infeasibility
    double tol_pr_, tol_du_;

    /// Minimum step size allowed
    double min_step_size_;

    /// Elastic mode
    bool elastic_mode_;

    /// Initial and maximum penalty parameter for elastic mode
    double gamma_0_, gamma_max_, gamma_1_min_;

    /// Initialize feasible qp's
    bool init_feasible_;

    /// Linesearch parameters
    ///@{
    double c1_;
    double beta_;
    casadi_int max_iter_ls_;
    casadi_int merit_memsize_;
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

    // Second order corrections
    bool so_corr_;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Generate code for the declarations of the C function */
    void codegen_declarations(CodeGenerator& g) const override;

    /// Access Conic
    const Function getConic() const { return qpsol_;}

    /// Print iteration header
    void print_iteration() const;

    /// Print iteration
    void print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf,
      double dx_norm, double rg, casadi_int ls_trials, bool ls_success,
      bool so_succes, std::string info) const;

    // Solve the QP subproblem: mode 0 = normal, mode 1 = SOC
    virtual int solve_QP(SqpmethodMemory* m, const double* H, const double* g,
      const double* lbdz, const double* ubdz,
      const double* A,
      double* x_opt, double* dlam, int mode) const;

    // Solve the QP subproblem for elastic mode
    virtual int solve_ela_QP(SqpmethodMemory* m, const double* H, const double* g,
      const double* lbdz, const double* ubdz, const double* A,
      double* x_opt, double* dlam) const;

    // Execute elastic mode: mode 0 = normal, mode 1 = SOC
    virtual int solve_elastic_mode(SqpmethodMemory* m, casadi_int* ela_it, double gamma_1,
      casadi_int ls_iter, bool ls_success, bool so_succes, double pr_inf,
      double du_inf, double dx_norminf, std::string* info, int mode) const;


    // Solve the QP subproblem
    void codegen_qp_solve(CodeGenerator& cg, const std::string& H, const std::string& g,
      const std::string& lbdz, const std::string& ubdz,
      const std::string& A, const std::string& x_opt, const std::string& dlam, int mode) const;

    // Solve the QP subproblem
    void codegen_qp_ela_solve(CodeGenerator& cg, const std::string& H, const std::string& g,
      const std::string& lbdz, const std::string& ubdz,
      const std::string& A, const std::string& x_opt, const std::string& dlam) const;

    // Execute elastic mode: mode 0 = normal, mode 1 = SOC
    void codegen_solve_elastic_mode(CodeGenerator& cg, int mode) const;

    // Codegen to calculate gama_1
    void codegen_calc_gamma_1(CodeGenerator& cg) const;


    // Calculate gamma_1
    double calc_gamma_1(SqpmethodMemory* m) const;

    /// A documentation string
    static const std::string meta_doc;


    /** \brief Serialize an object without type information */
    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize into MX */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new Sqpmethod(s); }

  protected:
    /** \brief Deserializing constructor */
    explicit Sqpmethod(DeserializingStream& s);

  private:
    void set_sqpmethod_prob();
  };

} // namespace casadi
/// \endcond
#endif // CASADI_SQPMETHOD_HPP
