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

  /// Input arguments of an ODE/DAE function [daeIn]
  enum DAEInput {
    /// Differential state [x]
    DAE_X,
    /// Algebraic state [z]
    DAE_Z,
    /// Parameter [p]
    DAE_P,
    /// Explicit time dependence [t]
    DAE_T,
    /// Number of arguments.
    DAE_NUM_IN
  };

  /// Output arguments of an DAE function [daeOut]
  enum DAEOutput {
    /// Right hand side of the implicit ODE [ode]
    DAE_ODE,
    /// Right hand side of algebraic equations [alg]
    DAE_ALG,
    /// Right hand side of quadratures equations [quad]
    DAE_QUAD,
    /// Number of arguments.
    DAE_NUM_OUT
  };

  /// Input arguments of an ODE/DAE backward integration function [rdaeIn]
  enum RDAEInput {
    /// Backward differential state [rx]
    RDAE_RX,
    /// Backward algebraic state [rz]
    RDAE_RZ,
    /// Backward  parameter vector [rp]
    RDAE_RP,
    /// Forward differential state [x]
    RDAE_X,
    /// Forward algebraic state [z]
    RDAE_Z,
    /// Parameter vector [p]
    RDAE_P,
    /// Explicit time dependence [t]
    RDAE_T,
    /// Number of arguments.
    RDAE_NUM_IN
  };

  /// Output arguments of an ODE/DAE backward integration function [rdaeOut]
  enum RDAEOutput {
    /// Right hand side of ODE. [ode]
    RDAE_ODE,
    /// Right hand side of algebraic equations. [alg]
    RDAE_ALG,
    /// Right hand side of quadratures. [quad]
    RDAE_QUAD,
    /// Number of arguments.
    RDAE_NUM_OUT
  };

  /// Input arguments of an integrator [integratorIn]
  enum IntegratorInput {
    /// Differential state at the initial time [x0]
    INTEGRATOR_X0,
    /// Parameters [p]
    INTEGRATOR_P,
    /// Initial guess for the algebraic variable [z0]
    INTEGRATOR_Z0,
    /// Backward differential state at the final time [rx0]
    INTEGRATOR_RX0,
    /// Backward parameter vector [rp]
    INTEGRATOR_RP,
    /// Initial guess for the backwards algebraic variable [rz0]
    INTEGRATOR_RZ0,
    /// Number of input arguments of an integrator
    INTEGRATOR_NUM_IN
  };

  /// Output arguments of an integrator [integratorOut]
  enum IntegratorOutput {
    /// Differential state at the final time [xf]
    INTEGRATOR_XF,
    /// Quadrature state at the final time [qf]
    INTEGRATOR_QF,
    /// Algebraic variable at the final time [zf]
    INTEGRATOR_ZF,
    /// Backward differential state at the initial time [rxf]
    INTEGRATOR_RXF,
    /// Backward quadrature state at the initial time [rqf]
    INTEGRATOR_RQF,
    /// Backward algebraic variable at the initial time [rzf]
    INTEGRATOR_RZF,
    /// Number of output arguments of an integrator
    INTEGRATOR_NUM_OUT
  };

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
  class CASADI_CORE_EXPORT Integrator : public Function {
  public:
    /// Default constructor
    Integrator();

    /** \brief  Integrator factory
    *
    * \param name \pluginargument{Integrator}
    * \param f dynamical system
    * \parblock
    * \copydoc scheme_DAEInput
    * \copydoc scheme_DAEOutput
    * \endparblock
    * \param g backwards system
    * \parblock
    * \copydoc scheme_RDAEInput
    * \copydoc scheme_RDAEOutput
    * \endparblock
    */
    Integrator(const std::string& name, const Function& f, const Function& g=Function());

    /// Clone
    Integrator clone() const;

    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;

    /// Access functions of the node
    IntegratorInternal* operator->();

    /// Access functions of the node
    const IntegratorInternal* operator->() const;

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

    /** \brief Generate a augmented DAE system with \a nfwd forward sensitivities
     *    and \a nadj adjoint sensitivities
     */
    std::pair<Function, Function> getAugmented(int nfwd, int nadj);

    /// Get the DAE
    Function getDAE();

    /// Set a stop time for the forward integration
    void setStopTime(double tf);

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };
} // namespace casadi

#endif // CASADI_INTEGRATOR_HPP
