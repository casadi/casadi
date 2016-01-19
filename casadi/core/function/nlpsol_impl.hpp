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


#ifndef CASADI_NLPSOL_IMPL_HPP
#define CASADI_NLPSOL_IMPL_HPP

#include "nlpsol.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL
namespace casadi {

  /** \brief Integrator memory */
  struct CASADI_EXPORT NlpsolMemory : public WorkMemory {
    // Outputs
    double *x, *f, *g, *lam_x, *lam_g, *lam_p;

    // Inputs
    const double *x0, *p, *lbx, *ubx, *lbg, *ubg, *lam_x0, *lam_g0;

    // Accumulated counts since last reset:
    int n_calc_f; // number of calls to calc_f
    int n_calc_g; // number of calls to calc_g
    int n_calc_grad_f; // number of calls to calc_grad_f
    int n_calc_jac_g; // number of calls to calc_jac_g
    int n_calc_hess_l; // number of calls to calc_hess_l
    int n_eval_callback; // number of calls to callback
    int n_iter; // number of iterations

    // Accumulated time since last reset:
    double t_calc_f; // time spent in calc_f
    double t_calc_g; // time spent in calc_g
    double t_calc_grad_f; // time spent in calc_grad_f
    double t_calc_jac_g; // time spent in calc_jac_g
    double t_calc_hess_l; // time spent in calc_hess_l

    /** \brief  Destructor */
    virtual ~NlpsolMemory() {}
  };

  /** \brief NLP solver storage class

      @copydoc Nlpsol_doc
      \author Joel Andersson
      \date 2010-2013
  */
  class CASADI_EXPORT
  Nlpsol : public FunctionInternal, public PluginInterface<Nlpsol> {

  public:
    /// Constructor
    Nlpsol(const std::string& name, const XProblem& nlp);

    /// Destructor
    virtual ~Nlpsol() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const { return NLPSOL_NUM_IN;}
    virtual size_t get_n_out() const { return NLPSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const;
    virtual Sparsity get_sparsity_out(int ind) const;
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::vector<std::string> get_ischeme() const { return nlpsol_in();}
    virtual std::vector<std::string> get_oscheme() const { return nlpsol_out();}
    /// @}

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual Memory* memory() const { return new NlpsolMemory();}

    /** \brief Initalize memory block */
    virtual void init_memory(Memory& mem) const;

    /** \brief Check if the inputs correspond to a well-posed problem */
    virtual void checkInputs(Memory& mem) const;

    /** \brief Get default input value */
    virtual double default_in(int ind) const;

    /// Number of variables
    int nx_;

    /// Number of constraints
    int ng_;

    /// Number of parameters
    int np_;

    /// callback function, executed at each iteration
    Function fcallback_;

    /// Execute the callback function only after this amount of iterations
    int callback_step_;

    // Evaluation errors are fatal
    bool eval_errors_fatal_;

    // Warn if initial bounds are violated
    bool warn_initial_bounds_;

    // Ignore errors in the iteration callbacks
    bool iteration_callback_ignore_errors_;

    /// The NLP
    XProblem nlp2_;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(Memory& mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    /** \brief Set the (temporary) work vectors */
    virtual void set_temp(Memory& mem, const double** arg, double** res,
                          int* iw, double* w) const;

    // Evaluate numerically
    virtual void eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const;

    // Solve the NLP
    virtual void solve(Memory& mem) const = 0;

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const XProblem& nlp);

    // No static functions exposed
    struct Exposed{ };

    // Calculate objective
    enum FIn { F_X, F_P, F_NUM_IN };
    enum FOut { F_F, F_NUM_OUT};
    Function f_fcn_;
    template<typename M> void _setup_f();
    void setup_f();
    int calc_f(NlpsolMemory& m, const double* x, const double* p, double* f) const;

    // Calculate constraints
    enum GIn { G_X, G_P, G_NUM_IN };
    enum GOut { G_G, G_NUM_OUT};
    Function g_fcn_;
    template<typename M> void _setup_g();
    void setup_g();
    int calc_g(NlpsolMemory& m, const double* x, const double* p, double* g) const;

    // Calculate both objective and constraints
    Function fg_fcn_;
    template<typename M> void _setup_fg();
    void setup_fg();
    int calc_fg(NlpsolMemory& m, const double* x, const double* p, double* f, double* g) const;

    // Calculate gradient of the objective
    enum GradFIn { GF_X, GF_P, GF_NUM_IN };
    enum GradFOut { GF_GF, GF_NUM_OUT};
    Function grad_f_fcn_;
    template<typename M> void _setup_grad_f();
    void setup_grad_f();
    int calc_grad_f(NlpsolMemory& m, const double* x,
                    const double* p, double* f, double* grad_f) const;

    // Calculate Jacobian of constraints
    enum JacGIn { JG_X, JG_P, JG_NUM_IN };
    enum JacGOut { JG_JG, JG_NUM_OUT};
    Function jac_g_fcn_;
    Sparsity jacg_sp_;
    template<typename M> void _setup_jac_g();
    void setup_jac_g();
    int calc_jac_g(NlpsolMemory& m, const double* x,
                   const double* p, double* g, double* jac_g) const;

    // Calculate Jacobian of gradient (note: sparse!)
    Function jac_f_fcn_;
    template<typename M> void _setup_jac_f();
    void setup_jac_f();
    int calc_jac_f(NlpsolMemory& m, const double* x,
                   const double* p, double* f, double* jac_f) const;

    // Calculate both gradient of the objective and Jacobian of constraints
    Function gf_jg_fcn_;
    template<typename M> void _setup_gf_jg();
    void setup_gf_jg();
    int calc_gf_jg(NlpsolMemory& m, const double* x,
                   const double* p, double* gf, double* jg) const;

    // Calculate Hessian of the Lagrangian constraints
    enum HessLagIn { HL_X, HL_P, HL_LAM_F, HL_LAM_G, HL_NUM_IN };
    enum HessLagOut { HL_HL, HL_NUM_OUT};
    Function hess_l_fcn_;
    Sparsity hesslag_sp_;
    template<typename M> void _setup_hess_l(bool tr, bool sym, bool diag);
    void setup_hess_l(bool tr=false, bool sym=false, bool diag=false);
    int calc_hess_l(NlpsolMemory& m, const double* x, const double* p,
                    const double* sigma, const double* lambda,
                    double* hl) const;

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "nlpsol";}

    // Get reduced Hessian
    virtual DM getReducedHessian();

    /// Read options from parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// WORKAROUND: Add an element to an std::vector stored in a GenericType:
    template<typename Type> static void append_to_vec(GenericType& t, Type el) {
      std::vector<Type> v = t;
      v.push_back(el);
      t = v;
    }

    /// Convert dictionary to Problem
    template<typename XType>
      static Problem<XType> map2problem(const std::map<std::string, XType>& d);

    /// Convert Problem to dictionary
    template<typename XType>
      static std::map<std::string, XType> problem2map(const Problem<XType>& d);

    /// Get the (legacy) dae forward function
    template<typename XType>
      static Problem<XType> fun2problem(Function nlp);
  };

  template<typename XType>
  Problem<XType> Nlpsol::map2problem(const std::map<std::string, XType>& d) {
    std::vector<XType> nl_in(NL_NUM_IN), nl_out(NL_NUM_OUT);
    for (auto&& i : d) {
      if (i.first=="x") {
        nl_in[NL_X]=i.second;
      } else if (i.first=="p") {
        nl_in[NL_P]=i.second;
      } else if (i.first=="f") {
        nl_out[NL_F]=i.second;
      } else if (i.first=="g") {
        nl_out[NL_G]=i.second;
      } else {
        casadi_error("No such field: " + i.first);
      }
    }
    return {nl_in, nl_out};
  }

  template<typename XType>
  std::map<std::string, XType> Nlpsol::problem2map(const Problem<XType>& d) {
    return {
        {"x", d.in[NL_X]},
        {"p", d.in[NL_P]},
        {"f", d.out[NL_F]},
        {"g", d.out[NL_G]},
      };
  }

  template<typename XType>
  Problem<XType> Nlpsol::fun2problem(Function nlp) {
    Problem<XType> p;
    p.in = XType::get_input(nlp);
    casadi_assert(p.in.size()==NL_NUM_IN);
    p.out = nlp(p.in);
    casadi_assert(p.out.size()==NL_NUM_OUT);
    return p;
  }

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_IMPL_HPP
