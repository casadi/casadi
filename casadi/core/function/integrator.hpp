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


#ifndef CASADI_INTEGRATOR_HPP
#define CASADI_INTEGRATOR_HPP

#include "function.hpp"
#include "linsol.hpp"
#include "rootfinder.hpp"

namespace casadi {

  /** \defgroup main_integrator
      Create an ODE/DAE integrator
      Solves an initial value problem (IVP) coupled to a terminal value problem
      with differential equation given as an implicit ODE coupled to an algebraic
      equation and a set of quadratures:

      \verbatim
      Initial conditions at t=t0
      x(t0)  = x0
      q(t0)  = 0

      Forward integration from t=t0 to t=tf
      der(x) = function(x, z, p, t)                  Forward ODE
      0 = fz(x, z, p, t)                  Forward algebraic equations
      der(q) = fq(x, z, p, t)                  Forward quadratures

      Terminal conditions at t=tf
      rx(tf)  = rx0
      rq(tf)  = 0

      Backward integration from t=tf to t=t0
      der(rx) = gx(rx, rz, rp, x, z, p, t)        Backward ODE
      0 = gz(rx, rz, rp, x, z, p, t)        Backward algebraic equations
      der(rq) = gq(rx, rz, rp, x, z, p, t)        Backward quadratures

      where we assume that both the forward and backwards integrations are index-1
      (i.e. dfz/dz, dgz/drz are invertible) and furthermore that
      gx, gz and gq have a linear dependency on rx, rz and rp.

      \endverbatim

      \generalsection{Integrator}
      \pluginssection{Integrator}

      \author Joel Andersson
      \date 2011-2015
  */
  /** \defgroup integrator
  * @copydoc main_integrator
  *  @{
  */

  /** \if EXPANDED
  * @copydoc main_integrator
  * \endif
  */
  ///@{
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const SXDict& dae, const Dict& opts=Dict());
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const MXDict& dae, const Dict& opts=Dict());
#ifndef SWIG
  CASADI_EXPORT Function integrator(const std::string& name, const std::string& solver,
                                    const Function& dae, const Dict& opts=Dict());
#endif // SWIG
  ///@}

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_integrator(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_integrator(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_integrator(const std::string& name);

  /** \brief Get input scheme of integrators */
  CASADI_EXPORT std::vector<std::string> integrator_in();

  /** \brief Get integrator output scheme of integrators */
  CASADI_EXPORT std::vector<std::string> integrator_out();

  /** \brief Get integrator input scheme name by index */
  CASADI_EXPORT std::string integrator_in(int ind);

  /** \brief Get output scheme name by index */
  CASADI_EXPORT std::string integrator_out(int ind);

  /** \brief Get the number of integrator inputs */
  CASADI_EXPORT int integrator_n_in();

  /** \brief Get the number of integrator outputs */
  CASADI_EXPORT int integrator_n_out();
  /** @} */

} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
