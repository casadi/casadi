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

#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "generic_integrator.hpp"

/** \defgroup DAE_doc
    Solves an initial value problem (IVP) coupled to a terminal value problem
    with differential equation given as an implicit ODE coupled to an algebraic
    equation and a set of quadratures:
    \verbatim
    Initial conditions at t=t0
    x(t0)  = x0
    q(t0)  = 0
  
    Forward integration from t=t0 to t=tf
    der(x) = fx(x,z,p,t)                  Forward ODE
    0 = fz(x,z,p,t)                  Forward algebraic equations
    der(q) = fq(x,z,p,t)                  Forward quadratures
  
    Terminal conditions at t=tf
    rx(tf)  = rx0
    rq(tf)  = 0
  
    Backward integration from t=tf to t=t0
    der(rx) = gx(rx,rz,rp,x,z,p,t)        Backward ODE
    0 = gz(rx,rz,rp,x,z,p,t)        Backward algebraic equations
    der(rq) = gq(rx,rz,rp,x,z,p,t)        Backward quadratures

    where we assume that both the forward and backwards integrations are index-1
    (i.e. dfz/dz, dgz/drz are invertible) and furthermore that 
    gx, gz and gq have a linear dependency on rx, rz and rp.
  
    \endverbatim 
*/
namespace CasADi{

  /// Input arguments of an ODE/DAE function [daeIn]
  enum DAEInput{
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
  enum DAEOutput{
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
  enum RDAEInput{
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
  enum RDAEOutput{
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
  enum IntegratorInput{
    /// Differential state at the initial time [x0]
    INTEGRATOR_X0, 
    /// Parameters [p]
    INTEGRATOR_P,
    /// Backward differential state at the final time [rx0]
    INTEGRATOR_RX0, 
    /// Backward parameter vector [rp]
    INTEGRATOR_RP, 
    /// Number of input arguments of an integrator
    INTEGRATOR_NUM_IN
  };

  /// Output arguments of an integrator [integratorOut]
  enum IntegratorOutput{
    ///  Differential state at the final time [xf]
    INTEGRATOR_XF,
    ///  Quadrature state at the final time [qf]
    INTEGRATOR_QF,
    ///  Backward differential state at the initial time [rxf]
    INTEGRATOR_RXF,
    ///  Backward quadrature state at the initial time [rqf]
    INTEGRATOR_RQF,
    /// Number of output arguments of an integrator
    INTEGRATOR_NUM_OUT
  };

  /// Forward declaration of internal class
  class IntegratorInternal;

  /** Integrator abstract base class
      @copydoc DAE_doc
  
      The Integrator class provides some additional functionality, such as getting the value of the state 
      and/or sensitivities at certain time points.
    
      The class does not specify the method used for the integration. This is defined in derived classes.

      \author Joel Andersson
      \date 2010
  */

  // grep "addOption" integrator_internal.cpp | perl -pe 's/addOption\((.*?),(.*?),(.*?)\);(.*\/\/ (.*))?/* \1 \2 \3 ...  \5\\n/'

  class Integrator : public GenericIntegrator{
  public:
    /// Default constructor
    Integrator();

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
     * \param nsens        Number of sensitivities to be propagated along with the integration forward in time
     * \param nsensB       Number of sensitivities to be propagated along with the integration backward in time
     * \param nsensB_store Number of sensitivities to be propagated along with the integration backward in time 
     *                     that depend on sensitivities propagated along with the integration forward in time
     */
    void reset(int nsens=0, int nsensB=0, int nsensB_store=0);

    /// Integrate forward until a specified time point 
    void integrate(double t_out);

    /** \brief Reset the backward problem
     * Time will be set to tf and backward state to input(INTEGRATOR_RX0)
     */
    void resetB();

    /// Integrate backward until a specified time point
    void integrateB(double t_out);

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /** \brief Generate a augmented DAE system with nfwd forward sensitivities and nadj adjoint sensitivities
     */
    std::pair<FX,FX> getAugmented(int nfwd, int nadj);
  
    /// Get the DAE
    FX getDAE();
  };

} // namespace CasADi

#endif //INTEGRATOR_HPP
