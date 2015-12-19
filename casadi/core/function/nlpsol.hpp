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


#ifndef CASADI_NLPSOL_HPP
#define CASADI_NLPSOL_HPP

#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

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

    /// Initialize
    virtual void init();

    /** \brief Check if the inputs correspond to a well-posed problem */
    virtual void checkInputs(void* mem) const;

    /// Set options that make the NLP solver more suitable for solving QPs
    virtual void setDefaultOptions(const std::string& recipe) {}

    /// Get or generate a function to calculate the gradient of the objective function
    virtual Function getGradF();

    /// Get or generate a function to calculate the jacobian of the objective function
    virtual Function getJacF();

    /// Get or generate a function to calculate the Jacobian of the constraint function
    virtual Function getJacG();

    /// Get or generate a function to calculate the gradient of the Lagrangian function
    virtual Function getGradLag();

    /// Get or generate a function to calculate the Hessian of the Lagrangian function
    virtual Function getHessLag();

    /// Get or generate the sparsity pattern of the Hessian of the Lagrangian
    virtual Sparsity getSpHessLag();

    /** \brief Get default input value */
    virtual double default_in(int ind) const;

    // Access the objective gradient function
    Function& gradF();

    // Access the objective jacobian function (sparse)
    Function& jacF();

    /// Access the Jacobian of the constraint function
    Function& jacG();

    /// Access the Hessian of the Lagrangian function
    Function& hessLag();

    /// Access the gradient of the Lagrangian function
    Function& gradLag();

    /// Get the sparsity pattern of the Hessian of the Lagrangian
    Sparsity& spHessLag();

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

    /// The NLP
    XProblem nlp2_;
    Function nlp_;

    // Gradient of the objective
    Function gradF_;

    // Gradient of the objective
    Function jacF_;

    // Jacobian of the constraints
    Function jacG_;

    // Hessian of the Lagrangian
    Function hessLag_;

    // Gradient of the Lagrangian
    Function gradLag_;

    // Sparsity pattern of the Hessian of the Lagrangian
    Sparsity spHessLag_;

    /// A reference to this object to be passed to the user functions
    Function ref_;

    // Inputs
    const double *x0_, *p_, *lbx_, *ubx_, *lbg_, *ubg_, *lam_x0_, *lam_g0_;

    // Get an element of the inputs
    inline double x0(int i) const { return x0_ ? x0_[i] : 0;}
    inline double p(int i) const { return p_ ? p_[i] : 0;}
    inline double lbx(int i) const { return lbx_ ? lbx_[i] : 0;}
    inline double ubx(int i) const { return ubx_ ? ubx_[i] : 0;}
    inline double lbg(int i) const { return lbg_ ? lbg_[i] : 0;}
    inline double ubg(int i) const { return ubg_ ? ubg_[i] : 0;}
    inline double lam_x0(int i) const { return lam_x0_ ? lam_x0_[i] : 0;}
    inline double lam_g0(int i) const { return lam_g0_ ? lam_g0_[i] : 0;}

    // Outputs
    double *x_, *f_, *g_, *lam_x_, *lam_g_, *lam_p_;

    // Work vectors
    const double** arg_;
    double** res_;
    int* iw_;
    double* w_;

    // Evaluate numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    // Reset the solver
    virtual void reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w);

    // Solve the NLP
    virtual void solve(void* mem) {}

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const XProblem& nlp);

    // No static functions exposed
    struct Exposed{ };

    // Calculate objective
    enum FIn { F_X, F_P, F_NUM_IN };
    enum FOut { F_F, F_NUM_OUT};
    Function f_fcn_;
    template<typename M> void setup_f();
    int calc_f(const double* x, const double* p, double* f);

    // Calculate constraints
    enum GIn { G_X, G_P, G_NUM_IN };
    enum GOut { G_G, G_NUM_OUT};
    Function g_fcn_;
    template<typename M> void setup_g();
    int calc_g(const double* x, const double* p, double* g);

    // Calculate both objective and constraints
    Function fg_fcn_;
    template<typename M> void setup_fg();
    int calc_fg(const double* x, const double* p, double* f, double* g);

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
    Sparsity jacg_sp_;
    template<typename M> void setup_jac_g();
    int calc_jac_g(double* jac_g=0);

    // Calculate both gradient of the objective and Jacobian of constraints
    Function gf_jg_fcn_;
    template<typename M> void setup_gf_jg();
    int calc_gf_jg(const double* x, const double* p, double* gf, double* jg);

    // Calculate Hessian of the Lagrangian constraints
    enum HessLagIn { HL_X, HL_P, HL_LAM_F, HL_LAM_G, HL_NUM_IN };
    enum HessLagOut { HL_HL, HL_NUM_OUT};
    Function hess_l_fcn_;
    Sparsity hesslag_sp_;
    template<typename M> void setup_hess_l();
    int calc_hess_l(const double* x, const double* p,
                    const double* sigma, const double* lambda,
                    double* hl);

    // Current solution
    double *xk_, lam_fk_, *lam_gk_, *lam_xk_;

    // Trigger recalculation?
    bool new_x_, new_lam_f_, new_lam_g_;

    // Current calculated quantities
    double fk_, *gk_, *grad_fk_, *jac_gk_, *hess_lk_, *grad_lk_;

    // Set primal variable
    void set_x(const double *x);

    // Set dual variable
    void set_lam_f(double lam_f);
    void set_lam_g(const double *lam_g);

    // Accumulated counts since last reset:
    int n_calc_f_; // number of calls to calc_f
    int n_calc_g_; // number of calls to calc_g
    int n_calc_grad_f_; // number of calls to calc_grad_f
    int n_calc_jac_g_; // number of calls to calc_jac_g
    int n_calc_hess_l_; // number of calls to calc_hess_l
    int n_eval_callback_; // number of calls to callback
    int n_iter_; // number of iterations

    // Accumulated time since last reset:
    double t_calc_f_; // time spent in calc_f
    double t_calc_g_; // time spent in calc_g
    double t_calc_grad_f_; // time spent in calc_grad_f
    double t_calc_jac_g_; // time spent in calc_jac_g
    double t_calc_hess_l_; // time spent in calc_hess_l

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
      static Function problem2fun(const Problem<XType>& d);

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
  Function Nlpsol::problem2fun(const Problem<XType>& d) {
    return Function("nlp", d.in, d.out, {"x", "p"}, {"f", "g"});
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

  template<typename M>
  void Nlpsol::setup_f() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg(F_NUM_IN);
    arg[F_X] = nlp.in[NL_X];
    arg[F_P] = nlp.in[NL_P];
    std::vector<M> res(F_NUM_OUT);
    res[F_F] = nlp.out[NL_F];
    f_fcn_ = Function("nlp_f", arg, res);
    alloc(f_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_g() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg(G_NUM_IN);
    arg[G_X] = nlp.in[NL_X];
    arg[G_P] = nlp.in[NL_P];
    std::vector<M> res(G_NUM_OUT);
    res[G_G] = nlp.out[NL_G];
    g_fcn_ = Function("nlp_g", arg, res);
    alloc(g_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_fg() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg = {nlp.in[NL_X], nlp.in[NL_P]};
    std::vector<M> res = {nlp.out[NL_F], nlp.out[NL_G]};
    fg_fcn_ = Function("nlp_fg", arg, res);
    alloc(fg_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_gf_jg() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg = {nlp.in[NL_X], nlp.in[NL_P]};
    std::vector<M> res = {M::gradient(nlp.out[NL_F], nlp.in[NL_X]),
                          M::jacobian(nlp.out[NL_G], nlp.in[NL_X])};
    gf_jg_fcn_ = Function("nlp_gf_jg", arg, res);
    alloc(gf_jg_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_grad_f() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg(GF_NUM_IN);
    arg[GF_X] = nlp.in[NL_X];
    arg[GF_P] = nlp.in[NL_P];
    std::vector<M> res(GF_NUM_OUT);
    res[GF_GF] = M::gradient(nlp.out[NL_F], nlp.in[NL_X]);
    grad_f_fcn_ = Function("nlp_grad_f", arg, res);
    alloc(grad_f_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_jac_g() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg(JG_NUM_IN);
    arg[JG_X] = nlp.in[NL_X];
    arg[JG_P] = nlp.in[NL_P];
    std::vector<M> res(JG_NUM_OUT);
    res[JG_JG] = M::jacobian(nlp.out[NL_G], nlp.in[NL_X]);
    jac_g_fcn_ = Function("nlp_jac_g", arg, res);
    jacg_sp_ = res[JG_JG].sparsity();
    alloc(jac_g_fcn_);
  }

  template<typename M>
  void Nlpsol::setup_hess_l() {
    const Problem<M>& nlp = nlp2_;
    std::vector<M> arg(HL_NUM_IN);
    M x = arg[HL_X] = nlp.in[NL_X];
    arg[HL_P] = nlp.in[NL_P];
    M f = nlp.out[NL_F];
    M g = nlp.out[NL_G];
    M lam_f = arg[HL_LAM_F] = M::sym("lam_f", f.sparsity());
    M lam_g = arg[HL_LAM_G] = M::sym("lam_g", g.sparsity());
    std::vector<M> res(HL_NUM_OUT);
    res[HL_HL] = triu(M::hessian(dot(lam_f, f) + dot(lam_g, g), x));
    hess_l_fcn_ = Function("nlp_hess_l", arg, res);
    hesslag_sp_ = res[HL_HL].sparsity();
    alloc(hess_l_fcn_);
  }

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_HPP
