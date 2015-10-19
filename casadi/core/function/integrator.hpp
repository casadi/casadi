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
#include "linear_solver.hpp"

/** \brief Base class for integrators
 *
 * \defgroup DAE_doc
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
*/
namespace casadi {
  /// Forward declaration of internal class
  class IntegratorInternal;

  // grep "addOption" integrator_internal.cpp |
  //   perl -pe 's/addOption\((.*?),(.*?), (.*?)\);(.*\/\/ (.*))?/* \1 \2 \3 ...  \5\\n/'

  /** Integrator abstract base class

      @copydoc DAE_doc

      The Integrator class provides some additional functionality, such as getting the value
      of the state and/or sensitivities at certain time points.

      \generalsection{Integrator}
      \pluginssection{Integrator}

      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT Integrator : public Function {
  public:
    /// Default constructor
    Integrator();

    ///@}
    /** \brief  Integrator factory (new syntax, includes initialization)
    *
    * \param solver \pluginargument{Integrator}
    * \param f dynamical system
    * \parblock
    * \copydoc scheme_DAEInput
    * \copydoc scheme_DAEOutput
    * \endparblock
    */
    Integrator(const std::string& name, const std::string& solver,
               const Function& f, const Dict& opts=Dict());
    Integrator(const std::string& name, const std::string& solver,
               const std::pair<Function, Function>& fg, const Dict& opts=Dict());
    ///@}

    /// Print solver statistics
    void printStats(std::ostream &stream=casadi::userOut()) const;

    /// Access functions of the node
    IntegratorInternal* operator->();

    /// Access functions of the node
    const IntegratorInternal* operator->() const;

#ifndef SWIG
    /** \brief Reset the forward problem
     * Time will be set to t0 and state to input(INTEGRATOR_X0)
     */
    void reset();

    /// Integrate forward until a specified time point
    void integrate(double t_out);

    /** \brief Reset the backward problem
     *
     * Time will be set to tf and backward state to input(INTEGRATOR_RX0)
     */
    void resetB();

    /// Integrate backward until a specified time point
    void integrateB(double t_out);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Set a stop time for the forward integration
    void setStopTime(double tf);
#endif // SWIG

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectNode* ptr);
  };
} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
