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

#ifndef CASADI_KNITRO_INTERNAL_HPP
#define CASADI_KNITRO_INTERNAL_HPP

#include <casadi/interfaces/knitro/casadi_nlpsolver_knitro_export.h>
#include <knitro.h>
#include "casadi/core/function/nlp_solver_internal.hpp"

/** \defgroup plugin_NlpSolver_knitro
  KNITRO interface
*/

/** \pluginsection{NlpSolver,knitro} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{NlpSolver,knitro}
     @copydoc NlpSolver_doc
     @copydoc plugin_NlpSolver_knitro
  */
  class CASADI_NLPSOLVER_KNITRO_EXPORT KnitroInternal : public NlpSolverInternal {

  public:
    explicit KnitroInternal(const Function& nlp);
    virtual ~KnitroInternal();
    virtual KnitroInternal* clone() const { return new KnitroInternal(*this);}

    /** \brief  Create a new NLP Solver */
    static NlpSolverInternal* creator(const Function& nlp)
    { return new KnitroInternal(nlp);}

    virtual void init();
    virtual void evaluate();

    // KNITRO callback functions
    void evalfc(const double* x, double& obj, double *c);
    void evalga(const double* x, double* objGrad, double* jac);
    void evalh(const double* x, const double* lambda, double* hessian);

    // KNITRO callback wrapper
    static int callback(const int evalRequestCode, const int n, const int m, const int nnzJ,
                        const int nnzH, const double * const x, const double * const lambda,
                        double * const obj, double * const c, double * const objGrad,
                        double * const jac, double * const hessian,
                        double * const hessVector, void *userParams);

    // KNITRO context pointer
    KTR_context_ptr kc_handle_;

    // KNITRO double parameter
    std::map<std::string, double> double_param_;

    // KNITRO int parameter
    std::map<std::string, int> int_param_;

    // KNITRO string parameter
    std::map<std::string, std::string> string_param_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_KNITRO_INTERNAL_HPP
