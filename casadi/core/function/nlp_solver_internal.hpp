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


#ifndef CASADI_NLP_SOLVER_INTERNAL_HPP
#define CASADI_NLP_SOLVER_INTERNAL_HPP

#include "nlp_solver.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"


/// \cond INTERNAL
namespace casadi {

/** \brief NLP solver storage class

  @copydoc NlpSolver_doc
  \author Joel Andersson
  \date 2010-2013
*/
  class CASADI_CORE_EXPORT
  NlpSolverInternal : public FunctionInternal,
                      public PluginInterface<NlpSolverInternal> {

  public:
    /// Constructor
    NlpSolverInternal(const Function& nlp);

    /// Destructor
    virtual ~NlpSolverInternal() = 0;

    /// Initialize
    virtual void init();

    /// Prints out a human readable report about possible constraint violations - all constraints
    void reportConstraints(std::ostream &stream=std::cout);

    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const;

    /// Warns the user about initial bounds, if option 'warn_initial_bounds' is true
    virtual void checkInitialBounds();

    /// Set options that make the NLP solver more suitable for solving QPs
    virtual void setQPOptions() {}

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
    Callback callback_;

    /// Execute the callback function only after this amount of iterations
    int callback_step_;

    // Evaluation errors are fatal
    bool eval_errors_fatal_;

    /// The NLP
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
    typedef NlpSolverInternal* (*Creator)(const Function& nlp);

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

  };

} // namespace casadi
/// \endcond
#endif // CASADI_NLP_SOLVER_INTERNAL_HPP
