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


#ifndef CASADI_KNITRO_INTERFACE_HPP
#define CASADI_KNITRO_INTERFACE_HPP

#include <casadi/interfaces/knitro/casadi_nlpsol_knitro_export.h>
#include <knitro.h>
#include "casadi/core/function/nlpsol.hpp"

/** \defgroup plugin_Nlpsol_knitro
  KNITRO interface
*/

/** \pluginsection{Nlpsol,knitro} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Nlpsol,knitro}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_knitro
  */
  class CASADI_NLPSOL_KNITRO_EXPORT KnitroInterface : public Nlpsol {

  public:
    explicit KnitroInterface(const std::string& name, const XProblem& nlp);
    virtual ~KnitroInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "knitro";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const XProblem& nlp) {
      return new KnitroInterface(name, nlp);
    }

    // Initialize the solver
    virtual void init();

    // Reset the solver
    virtual void reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w);

    // Solve the NLP
    virtual void solve(void* mem);

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

    // Inputs
    double *wx_, *wlbx_, *wubx_, *wlbg_, *wubg_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_KNITRO_INTERFACE_HPP
