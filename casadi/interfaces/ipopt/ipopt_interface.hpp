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


#ifndef CASADI_IPOPT_INTERFACE_HPP
#define CASADI_IPOPT_INTERFACE_HPP

#include <casadi/interfaces/ipopt/casadi_nlpsol_ipopt_export.h>
#include "casadi/core/function/nlpsol.hpp"
#include "casadi/core/timing.hpp"

/** \defgroup plugin_Nlpsol_ipopt
 *
 * When in warmstart mode, output NLPSOL_LAM_X may be used as input
 *
 * NOTE: Even when max_iter == 0, it is not guaranteed that
 * input(NLPSOL_X0) == output(NLPSOL_X).
 * Indeed if bounds on X or constraints are unmet, they will differ.
 *
 * For a good tutorial on IPOPT, see
 * http://drops.dagstuhl.de/volltexte/2009/2089/pdf/09061.WaechterAndreas.Paper.2089.pdf
 *
 * A good resource about the algorithms in IPOPT is: Wachter and L. T. Biegler,
 * On the Implementation of an Interior-Point Filter Line-Search Algorithm for
 * Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57,
 * 2006 (As Research Report RC 23149, IBM T. J. Watson Research Center, Yorktown, USA
 *
 * Caveats:
 * * with default options, multipliers for the decision variables are wrong for equality
 * constraints.
 * Change the 'fixed_variable_treatment' to 'make_constraint' or 'relax_bounds' to obtain
 * correct results.
 *
 */

/** \pluginsection{Nlpsol,ipopt} **/

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Nlpsol,ipopt}

      @copydoc Nlpsol_doc
      @copydoc plugin_Nlpsol_ipopt
  */
  class CASADI_NLPSOL_IPOPT_EXPORT IpoptInterface : public Nlpsol {
    friend class IpoptUserClass;

  public:
    explicit IpoptInterface(const std::string& name, const XProblem& nlp);
    virtual ~IpoptInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "ipopt";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const XProblem& nlp) {
      return new IpoptInterface(name, nlp);
    }

    // Initialize the solver
    virtual void init();

    // Reset the solver
    virtual void reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w);

    // Solve the NLP
    virtual void solve(void* mem);

    /// Set default options for a given recipe
    virtual void setDefaultOptions(const std::vector<std::string>& recipes);

    // Get reduced Hessian
    virtual DMatrix getReducedHessian();

    /// Exact Hessian?
    bool exact_hessian_;

    /** NOTE:
     * To allow this header file to be free of IPOPT types
     * (that are sometimes declared outside their scope!) and after
     * experiencing problems with working with IPOPT classes without
     * IPOPT smart pointers, we work with dynamically allocated IPOPT
     * smart pointers in this interface, that are stored as void
     * pointers in the interface.
     *
     */
    void *userclass_;
    void* app_;
#ifdef WITH_SIPOPT
    void* app_sens_;
#endif // WITH_SIPOPT

    /// All IPOPT options
    std::map<std::string, TypeID> ops_;

    // Sparsities
    Sparsity jacg_sp_, hesslag_sp_;

    // Current solution
    double *xk_, lam_fk_, *lam_gk_, *lam_xk_;

    // Trigger recalculation?
    bool new_x_, new_lam_f_, new_lam_g_;

    // Current calculated quantities
    double fk_, *gk_, *grad_fk_, *jac_gk_, *hess_lk_, *grad_lk_;

    // What is up-to-date?
    bool have_fk_, have_gk_, have_grad_fk_, have_jac_gk_, have_hess_lk_, have_grad_lk_;

    // Set primal variable
    void set_x(const double *x);

    // Set dual variable
    void set_lam_f(double lam_f);
    void set_lam_g(const double *lam_g);

    // Calculate objective
    enum FIn { F_X, F_P, F_NUM_IN };
    enum FOut { F_F, F_NUM_OUT};
    Function f_fcn_;
    template<typename M> void setup_f();
    int calc_f(double* f=0);

    // Calculate constraints
    enum GIn { G_X, G_P, G_NUM_IN };
    enum GOut { G_G, G_NUM_OUT};
    Function g_fcn_;
    template<typename M> void setup_g();
    int calc_g(double* g=0);

    // Calculate gradient of the objective
    enum GradFIn { GF_X, GF_P, GF_NUM_IN };
    enum GradFOut { GF_GF, GF_NUM_OUT};
    Function grad_f_fcn_;
    template<typename M> void setup_grad_f();
    int calc_grad_f(double* grad_f=0);

    // Calculate Jacobian of constraints
    enum JacGIn { JG_X, JG_P, JG_NUM_IN };
    enum JacGOut { JG_JG, JG_NUM_OUT};
    Function jac_g_fcn_;
    template<typename M> void setup_jac_g();
    int calc_jac_g(double* jac_g=0);

    // Calculate Hessian of the Lagrangian constraints
    enum HessLagIn { HL_X, HL_P, HL_LAM_F, HL_LAM_G, HL_NUM_IN };
    enum HessLagOut { HL_HL, HL_NUM_OUT};
    Function hess_l_fcn_;
    template<typename M> void setup_hess_l();
    int calc_hess_l(double* hess_l=0);

    // Setup all functions
    template<typename M> void setup();

    // Ipopt callback functions
    void finalize_solution(const double* x, const double* z_L, const double* z_U, const double* g,
                           const double* lambda, double obj_value, int iter_count);
    bool get_bounds_info(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u);
    bool get_starting_point(int n, bool init_x, double* x, bool init_z, double* z_L, double* z_U,
                            int m, bool init_lambda, double* lambda);
    void get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag);
    int get_number_of_nonlinear_variables();
    bool get_list_of_nonlinear_variables(int num_nonlin_vars, int* pos_nonlin_vars);
    bool intermediate_callback(const double* x, const double* z_L, const double* z_U,
                               const double* g,
                               const double* lambda, double obj_value, int iter,
                               double inf_pr, double inf_du, double mu, double d_norm,
                               double regularization_size, double alpha_du, double alpha_pr,
                               int ls_trials, bool full_callback);
    bool get_var_con_metadata(int n,
                              std::map<std::string, std::vector<std::string> >& var_string_md,
                              std::map<std::string, std::vector<int> >& var_integer_md,
                              std::map<std::string, std::vector<double> >& var_numeric_md,
                              int m,
                              std::map<std::string, std::vector<std::string> >& con_string_md,
                              std::map<std::string, std::vector<int> >& con_integer_md,
                              std::map<std::string, std::vector<double> >& con_numeric_md);

    void finalize_metadata(int n,
                           const std::map<std::string, std::vector<std::string> >& var_string_md,
                           const std::map<std::string, std::vector<int> >& var_integer_md,
                           const std::map<std::string, std::vector<double> >& var_numeric_md,
                           int m,
                           const std::map<std::string, std::vector<std::string> >& con_string_md,
                           const std::map<std::string, std::vector<int> >& con_integer_md,
                           const std::map<std::string, std::vector<double> >& con_numeric_md);

    // Accumulated time since last reset:
    double t_calc_f_; // time spent in calc_f
    double t_calc_g_; // time spent in calc_g
    double t_calc_grad_f_; // time spent in calc_grad_f
    double t_calc_jac_g_; // time spent in calc_jac_g
    double t_calc_hess_l_; // time spent in calc_hess_l

    // Timings for different parts of the main loop
    DiffTime t_callback_fun_;  // time spent in callback function
    DiffTime t_callback_prepare_; // time spent in callback preparation
    DiffTime t_mainloop_; // time spent in the main loop of the solver

    // Accumulated counts since last reset:
    int n_calc_f_; // number of calls to calc_f
    int n_calc_g_; // number of calls to calc_g
    int n_calc_grad_f_; // number of calls to calc_grad_f
    int n_calc_jac_g_; // number of calls to calc_jac_g
    int n_calc_hess_l_; // number of calls to calc_hess_l
    int n_eval_callback_; // number of calls to callback
    int n_iter_; // number of iterations

    // For parametric sensitivities with sIPOPT
#ifdef WITH_SIPOPT
    bool run_sens_;
    bool compute_red_hessian_;
    DMatrix red_hess_;
#endif // WITH_SIPOPT

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_IPOPT_INTERFACE_HPP
