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


#ifndef CASADI_WORHP_INTERFACE_HPP
#define CASADI_WORHP_INTERFACE_HPP

#include "casadi/core/function/nlp_solver_internal.hpp"
#include <casadi/interfaces/worhp/casadi_nlpsolver_worhp_export.h>

// GCC_VERSION is defined in 'worhp.h'
#ifdef GCC_VERSION
#undef GCC_VERSION
#endif

// Workaround for Clang, but should not be a problem for other compilers, #771
#define _Bool bool

#include <worhp.h>

// MACROs that pollute our code
#undef Q
/**\defgroup plugin_NlpSolver_worhp
 WORHP interface
*/
/** \pluginsection{NlpSolver,worhp} **/

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{NlpSolver,worhp}
     @copydoc NlpSolver_doc
     @copydoc plugin_NlpSolver_worhp
  */
  class CASADI_NLPSOLVER_WORHP_EXPORT WorhpInterface : public NlpSolverInternal {

  public:
    // Constructor
    explicit WorhpInterface(const Function& nlp);

    // Destructor
    virtual ~WorhpInterface();

    // Clone function
    virtual WorhpInterface* clone() const;

    /** \brief  Create a new NLP Solver */
    static NlpSolverInternal* creator(const Function& nlp)
    { return new WorhpInterface(nlp);}

    // Reset solver
    void reset();

    // (Re)initialize
    virtual void init();

    // Solve the NLP
    virtual void evaluate();


    virtual void setQPOptions();

    /// Read options from worhp parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// Exact Hessian?
    bool exact_hessian_;

    OptVar    worhp_o_;
    Workspace worhp_w_;
    Params    worhp_p_;
    Control   worhp_c_;

    std::map<int, std::string> status_;
    std::map<std::string, opt_type> ops_;

    // Worhp callback functions
    bool eval_f(const double* x, double scale, double& obj_value);
    bool eval_grad_f(const double* x, double scale , double* grad_f);
    bool eval_g(const double* x, double* g);
    bool eval_jac_g(const double* x, double* values);
    bool eval_h(const double* x, double obj_factor, const double* lambda, double* values);

    // Accumulated time since last reset:
    double t_eval_f_; // time spent in eval_f
    double t_eval_grad_f_; // time spent in eval_grad_f
    double t_eval_g_; // time spent in eval_g
    double t_eval_jac_g_; // time spent in eval_jac_g
    double t_eval_h_; // time spent in eval_h
    double t_callback_fun_;  // time spent in callback function
    double t_callback_prepare_; // time spent in callback preparation
    double t_mainloop_; // time spent in the main loop of the solver

    // Accumulated counts since last reset:
    int n_eval_f_; // number of calls to eval_f
    int n_eval_grad_f_; // number of calls to eval_grad_f
    int n_eval_g_; // number of calls to eval_g
    int n_eval_jac_g_; // number of calls to eval_jac_g
    int n_eval_h_; // number of calls to eval_h

    std::string formatStatus(int status) const;

    /// Pass the supplied options to Worhp
    void passOptions();

    // Calculate the status message map
    static std::map<int, std::string> calc_flagmap();

    // Error status map
    static std::map<int, std::string> flagmap;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_WORHP_INTERFACE_HPP
