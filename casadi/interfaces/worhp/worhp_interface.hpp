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

#include "casadi/core/function/nlpsol.hpp"
#include <casadi/interfaces/worhp/casadi_nlpsol_worhp_export.h>

// GCC_VERSION is defined in 'worhp.h'
#ifdef GCC_VERSION
#undef GCC_VERSION
#endif

// Workaround for Clang, but should not be a problem for other compilers, #771
#define _Bool bool

#include <worhp.h>

// MACROs that pollute our code
#undef Q
/**\defgroup plugin_Nlpsol_worhp
 WORHP interface
*/
/** \pluginsection{Nlpsol,worhp} **/

/// \cond INTERNAL
namespace casadi {

  struct CASADI_NLPSOL_WORHP_EXPORT WorhpMemory : public NlpsolMemory {
    OptVar    worhp_o;
    Workspace worhp_w;
    Params    worhp_p;
    Control   worhp_c;

    // Accumulated time since last reset:
    double t_eval_f; // time spent in eval_f
    double t_eval_grad_f; // time spent in eval_grad_f
    double t_eval_g; // time spent in eval_g
    double t_eval_jac_g; // time spent in eval_jac_g
    double t_eval_h; // time spent in eval_h
    double t_callback_fun;  // time spent in callback function
    double t_callback_prepare; // time spent in callback preparation
    double t_mainloop; // time spent in the main loop of the solver

    // Accumulated counts since last reset:
    int n_eval_f; // number of calls to eval_f
    int n_eval_grad_f; // number of calls to eval_grad_f
    int n_eval_g; // number of calls to eval_g
    int n_eval_jac_g; // number of calls to eval_jac_g
    int n_eval_h; // number of calls to eval_h

    // Stats
    int iter;
    int iter_sqp;
    double inf_pr;
    double inf_du;
    double alpha_pr;
    int return_code;
    const char* return_status;

    /// Constructor
    WorhpMemory();

    /// Destructor
    virtual ~WorhpMemory();
  };

  /** \brief \pluginbrief{Nlpsol,worhp}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_worhp
  */
  class CASADI_NLPSOL_WORHP_EXPORT WorhpInterface : public Nlpsol {

  public:
    // Constructor
    explicit WorhpInterface(const std::string& name, const XProblem& nlp);

    // Destructor
    virtual ~WorhpInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "worhp";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const XProblem& nlp) {
      return new WorhpInterface(name, nlp);
    }

    // Reset solver
    void reset();

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual Memory* memory() const { return new WorhpMemory();}

    /** \brief Initalize memory block */
    virtual void init_memory(Memory& mem) const;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(Memory& mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Solve the NLP
    virtual void solve(Memory& mem) const;

    // Options
    std::map<std::string, bool> bool_opts_;
    std::map<std::string, int> int_opts_;
    std::map<std::string, double> double_opts_;
    bool print_time_;

    // WORHP return codes
    static const char* return_codes(int flag);

    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_WORHP_INTERFACE_HPP
