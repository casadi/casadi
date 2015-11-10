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
    virtual const double& default_in(int ind) const;

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
    virtual void evalD(const double** arg, double** res, int* iw, double* w, void* mem);

    // Reset the solver
    virtual void reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w);

    // Solve the NLP
    virtual void solve(void* mem) {}

    // Creator function for internal class
    typedef Nlpsol* (*Creator)(const std::string& name, const XProblem& nlp);

    // No static functions exposed
    struct Exposed{ };

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

} // namespace casadi
/// \endcond
#endif // CASADI_NLPSOL_HPP
