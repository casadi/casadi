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


#ifndef CASADI_RK_IVPSOL_HPP
#define CASADI_RK_IVPSOL_HPP

#include "fixed_step_integrator.hpp"
#include <casadi/solvers/casadi_ivpsol_rk_export.h>

/** \defgroup plugin_Ivpsol_rk
      Fixed-step explicit Runge-Kutta integrator for ODEs
      Currently implements RK4.

      The method is still under development
*/
/** \pluginsection{Ivpsol,rk} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Ivpsol,rk}


      @copydoc DAE_doc
      @copydoc plugin_Ivpsol_rk

      \author Joel Andersson
      \date 2011-2014
  */
  class CASADI_IVPSOL_RK_EXPORT RkIvpsol : public FixedStepIvpsol {
  public:

    /// Constructor
    explicit RkIvpsol(const std::string& name, const XProblem& dae);

    /** \brief  Create a new integrator */
    static Ivpsol* creator(const std::string& name, const XProblem& dae) {
      return new RkIvpsol(name, dae);
    }

    /// Destructor
    virtual ~RkIvpsol();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "rk";}

    /// Initialize stage
    virtual void init();

    /// Setup F and G
    virtual void setupFG();

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_RK_IVPSOL_HPP
