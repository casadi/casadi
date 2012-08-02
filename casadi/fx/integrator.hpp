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

#include "fx.hpp"
#include "linear_solver.hpp"
//#define NEW_INTEGRATOR

/** \defgroup DAE_doc
  Solves the following initial value problem (IVP):
  
  \verbatim
  F(t,x,der(x),z,p) == 0
  x(t0) = x0
  over a time interval [t0, tf].
  \endverbatim 
  
  NOTE: The ODE/DAE initial-value problem formulation in CasADi is being replaced with a new 
  formulation which solves an initial value problem (IVP) coupled to a terminal value problem
  with differential equation given as an implicit ODE coupled to an algebraic
  equation and a set of quadratures:
  \verbatim
  Initial conditions at t=t0
    x(t0)  = x0
    q(t0)  = 0
  
  Forward integration from t=t0 to t=tf
         0 = fx(x,z,p,t,der(x))           Forward ODE
         0 = fz(x,z,p,t)                  Forward algebraic equations
    der(q) = fq(x,z,p,t)                  Forward quadratures
  
  Terminal conditions at t=tf
    rx(tf)  = rx0
    rq(tf)  = 0
  
  Backward integration from t=tf to t=t0
          0 = gx(rx,rz,x,z,p,t,der(rx))   Backward ODE
          0 = gz(rx,rz,x,z,p,t)           Backward algebraic equations
    der(rq) = gq(rx,rz,x,z,p,t)           Backward quadratures

  where we assume that both the forward and backwards integrations are index-1
  (i.e. dfx/dxdot, dfz/dz, dgz/drz, dgx/drxdot are invertible) and furthermore that 
  gx, gz and gq have a linear dependency on rx and rz and that f_x and g_x have a 
  linear dependence on xdot and rxdot respectively.
  \endverbatim 
*/

/** \defgroup ODE_doc
  Solves the following initial value problem (IVP):
  
  \verbatim
  xdot == f(t,x,p)
  from t0 to tf
  
  given the initial condition
  x(t0) == x0;
  \endverbatim 
*/

namespace CasADi{

/// Input arguments of an ODE/DAE function
enum DAEInput{
  /** Differential state */
  DAE_X,
  /** Algebraic state */
  DAE_Z,
  /** Parameter */
  DAE_P,
  /** Explicit time dependence */
  DAE_T,
  /** Time derivative of differential states */
  DAE_XDOT,
  /** Number of arguments. */
  DAE_NUM_IN
};

/// Helper function to create ODE/DAE forward integration function input arguments
template<class M>
std::vector<M> daeIn(const M& x, const M& z=M(), const M& p=M(), const M& t=M(), const M& xdot=M()){
  M ret[DAE_NUM_IN] = {x,z,p,t,xdot};
  return std::vector<M>(ret,ret+DAE_NUM_IN);
}
#ifdef SWIG
%template(daeIn) daeIn<SXMatrix>;
%template(daeIn) daeIn<MX>;
#endif //SWIG

/// Output arguments of an DAE function
enum DAEOutput{
  /** Right hand side of the implicit ODE */
  DAE_ODE,
  /** Right hand side of algebraic equations */
  DAE_ALG,
  /** Right hand side of quadratures equations */
  DAE_QUAD,
  /** Number of arguments. */
  DAE_NUM_OUT
};

/// Helper function to create DAE forward integration function output arguments
template<class M>
std::vector<M> daeOut(const M& ode, const M& alg=M(), const M& quad=M()){
  M ret[DAE_NUM_OUT] = {ode,alg,quad};
  return std::vector<M>(ret,ret+DAE_NUM_OUT);
}
#ifdef SWIG
%template(daeOut) daeOut<SXMatrix>;
%template(daeOut) daeOut<MX>;
#endif //SWIG

/// Input arguments of an ODE/DAE backward integration function 
enum RDAEInput{
  /** Backward differential state */
  RDAE_RX,
  /** Backward algebraic state */
  RDAE_RZ,
  /** Forward differential state */
  RDAE_X,
  /** Forward algebraic state */
  RDAE_Z,
  /** Parameter vector */
  RDAE_P,
  /** Explicit time dependence */
  RDAE_T,
  /** Time derivative of differential states */
  RDAE_XDOT,
  /** Time derivative of backward differential state */
  RDAE_RXDOT,
  /** Number of arguments. */
  RDAE_NUM_IN
};

/// Helper function to create ODE/DAE backward integration function input arguments
template<class M>
std::vector<M> rdaeIn(const M& rx, const M& rz=M(), const M& x=M(), const M& z=M(), const M& p=M(), const M& t=M(), const M& xdot=M(), const M& rxdot=M()){
  M ret[RDAE_NUM_IN] = {rx,rz,x,z,p,t,xdot,rxdot};
  return std::vector<M>(ret,ret+RDAE_NUM_IN);
}
#ifdef SWIG
%template(rdaeIn) rdaeIn<SXMatrix>;
%template(rdaeIn) rdaeIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE backward integration function
enum RDAEOutput{
  /** Right hand side of ODE.*/
  RDAE_ODE,
  /** Right hand side of algebraic equations.*/
  RDAE_ALG,
  /** Right hand side of quadratures.*/
  RDAE_QUAD,
  /** Number of arguments. */
  RDAE_NUM_OUT
};

/// Helper function to create ODE/DAE backward integration function output arguments
template<class M>
std::vector<M> rdaeOut(const M& ode, const M& alg=M(), const M& quad=M()){
  M ret[RDAE_NUM_OUT] = {ode,alg,quad};
  return std::vector<M>(ret,ret+RDAE_NUM_OUT);
}
#ifdef SWIG
%template(rdaeOut) rdaeOut<SXMatrix>;
%template(rdaeOut) rdaeOut<MX>;
#endif //SWIG

/// Input arguments of an integrator
enum IntegratorInput{
  /** Differential state at the initial time */
  INTEGRATOR_X0, 
  /** Parameters */
  INTEGRATOR_P,
  /** Backward differential state at the final time */
  INTEGRATOR_RX0, 
  /** Number of input arguments of an integrator */
  INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum IntegratorOutput{
  /**  Differential state at the final time */
  INTEGRATOR_XF,
  /**  Quadrature state at the final time */
  INTEGRATOR_QF,
  /**  Backward differential state at the initial time */
  INTEGRATOR_RXF,
  /**  Backward quadrature state at the initial time */
  INTEGRATOR_RQF,
    /** Number of output arguments of an integrator */
  INTEGRATOR_NUM_OUT
};

/// Forward declaration of internal class
class IntegratorInternal;

/** Integrator abstract base class
  @copydoc DAE_doc
  
  The Integrator class provides some additional functionality, such as getting the value of the state and/or sensitivities at certain time points.
  Controls are assumed to be parametrized at this point. 
    
  The class does not specify how the function F above should be represented, nor the method used for the integration, but assumes that it steps forward in time 
  (ruling out collocation in particular). The actual form of the ODE/DAE is defined in the derived classes.

  \author Joel Andersson
  \date 2010
*/

// grep "addOption" integrator_internal.cpp | perl -pe 's/addOption\((.*?),(.*?),(.*?)\);(.*\/\/ (.*))?/* \1 \2 \3 ...  \5\\n/'

class Integrator : public FX{
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
  
  /// Reset the solver and bring the time back to t0 and state back to INTEGRATOR_X0
  void reset(int nfdir=0);

  /// Integrate until a specified time point 
  void integrate(double t_out);

  /// Reset the solver of the adjoint problem and take time to tf
  void resetAdj();

  /// Integrate backwards in time until a specified time point
  void integrateAdj(double t_out);

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
