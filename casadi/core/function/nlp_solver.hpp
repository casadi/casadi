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


#ifndef CASADI_NLP_SOLVER_HPP
#define CASADI_NLP_SOLVER_HPP

#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief NLP solver storage class

  @copydoc NlpSolver_doc
  \author Joel Andersson
  \date 2010-2013
*/
  class CASADI_EXPORT
  NlpSolverInternal : public FunctionInternal,
                      public PluginInterface<NlpSolverInternal> {

  public:
    /// Constructor
    NlpSolverInternal(const std::string& name, const XProblem& nlp);

    /// Destructor
    virtual ~NlpSolverInternal() = 0;

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() const { return NLP_SOLVER_NUM_IN;}
    virtual size_t get_n_out() const { return NLP_SOLVER_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int ind) const;
    virtual Sparsity get_sparsity_out(int ind) const;
    /// @}

    /// Initialize
    virtual void init();

    /// Prints out a human readable report about possible constraint violations - all constraints
    void reportConstraints(std::ostream &stream=casadi::userOut());

    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const;

    /// Warns the user about initial bounds, if option 'warn_initial_bounds' is true
    virtual void checkInitialBounds();

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

    // Creator function for internal class
    typedef NlpSolverInternal* (*Creator)(const std::string& name, const XProblem& nlp);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "nlp";}

    // Get reduced Hessian
    virtual DMatrix getReducedHessian();

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
  Problem<XType> NlpSolverInternal::map2problem(const std::map<std::string, XType>& d) {
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
  std::map<std::string, XType> NlpSolverInternal::problem2map(const Problem<XType>& d) {
    return {
        {"x", d.in[NL_X]},
        {"p", d.in[NL_P]},
        {"f", d.out[NL_F]},
        {"g", d.out[NL_G]},
      };
  }

  template<typename XType>
  Function NlpSolverInternal::problem2fun(const Problem<XType>& d) {
    return Function("nlp", d.in, d.out, {"x", "p"}, {"f", "g"});
  }

  template<typename XType>
  Problem<XType> NlpSolverInternal::fun2problem(Function nlp) {
    Problem<XType> p;
    p.in = XType::get_input(nlp);
    casadi_assert(p.in.size()==NL_NUM_IN);
    p.out = nlp(p.in);
    casadi_assert(p.out.size()==NL_NUM_OUT);
    return p;
  }

} // namespace casadi
/// \endcond
#endif // CASADI_NLP_SOLVER_HPP
