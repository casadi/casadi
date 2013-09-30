/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef WORHP_INTERNAL_HPP
#define WORHP_INTERNAL_HPP

#include "symbolic/fx/nlp_solver_internal.hpp"

// GCC_VERSION is defined in worhp.h
#ifdef GCC_VERSION
#undef GCC_VERSION
#endif

// Workaround for Clang, but should not be a problem for other compilers, #771
#define _Bool bool

#include "worhp.h"

namespace CasADi{

  /**
     @copydoc NLPSolver_doc
  */
  class WorhpInternal : public NLPSolverInternal{

  public:
    // Constructor
    explicit WorhpInternal(const FX& nlp);

    // Destructor
    virtual ~WorhpInternal();

    // Clone function
    virtual WorhpInternal* clone() const;

    // Reset solver
    void reset();

    // (Re)initialize
    virtual void init();
    
    // Solve the NLP
    virtual void evaluate(int nfdir, int nadir);

    
    virtual void setQPOptions();

    /// Read options from worhp parameter xml
    void setOptionsFromFile(const std::string & file);
  
    /// Exact Hessian?
    bool exact_hessian_;

    // Sparsity pattern of the transpose of jacG
    CRSSparsity spJacG_T_;
    std::vector<int> jacG_tmp_;

    OptVar    worhp_o_;
    Workspace worhp_w_;
    Params    worhp_p_;
    Control   worhp_c_;
  
    std::map<int,std::string> status_;
    std::map<std::string,opt_type> ops_;

    // Worhp callback functions
    bool eval_f(const double* x, double scale, double& obj_value);
    bool eval_grad_f(const double* x, double scale , double* grad_f);
    bool eval_g(const double* x, double* g);
    bool eval_jac_g(const double* x, double* values);
    bool eval_h(const double* x, double obj_factor, const double* lambda, double* values);
  
    // Accummulated time since last reset:
    double t_eval_f_; // time spent in eval_f
    double t_eval_grad_f_; // time spent in eval_grad_f
    double t_eval_g_; // time spent in eval_g
    double t_eval_jac_g_; // time spent in eval_jac_g
    double t_eval_h_; // time spent in eval_h
    double t_callback_fun_;  // time spent in callback function
    double t_callback_prepare_; // time spent in callback preparation
    double t_mainloop_; // time spent in the main loop of the solver
  
    std::string formatStatus(int status) const;
  
    /// Pass the supplied options to Worhp
    void passOptions();

    // Calculate the status message map
    static std::map<int,std::string> calc_flagmap();

    // Error status map
    static std::map<int,std::string> flagmap;
  
  };

} // namespace CasADi

#endif //WORHP_INTERNAL_HPP
