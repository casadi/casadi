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


/** \defgroup DAE_doc
  Solves the following initial value problem (IVP):
  
  \verbatim
  F(t,x,der(x),z,p) == 0
  x(t0) = x0
  over a time interval [t0, tf].
  \endverbatim 
  
  NOTE: The ODE/DAE formulation in CasADi will be replaced by a more general
  forward/backward two-point boundary value problem with the differential equation
  given as a semi-explicit DAE with quadrature states:
  \verbatim
  Initial conditions at t=t0
    x_d(t0)  = x_d0
    x_q(t0)  = x_q0
  
  Forward integration from t=t0 to t=tf
    der(x_d) = f_d(t,x_d,x_a,p)           Forward ODE
    der(x_q) = f_q(t,x_d,x_a,p)           Forward quadratures
    der(x_a) = f_a(t,x_d,x_a,p)           Forward algebraic equations
  
  Terminal conditions at t=tf
    y_d(tf)  = h_d(x_d(tf),x_q(tf),p)
    y_q(tf)  = h_q(x_d(tf),x_q(tf),p)
  
  Backward integration from t=tf to t=t0
    der(y_d) = g_d(t,x_d,x_a,y_d,y_a,p)   Backward ODE
    der(y_q) = g_q(t,x_d,x_a,y_d,y_a,p)   Backward quadratures
    der(y_a) = g_a(t,x_d,x_a,y_d,y_a,p)   Backward algebraic equations

  Where we assume that both the forward and backwards integrations are index-1
  (i.e. d_f_a/d_x_a and d_q_a/d_y_a are invertible) and furthermore that 
  g_d, g_q and g_a have a linear dependency on y_d and y_a.
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
  /** Time. (1-by-1) */
  DAE_T,
  /** State vector (matrix). Should have same amount of non-zeros as DAEOutput:DAE_RES */
  DAE_Y,
  /** Parameter vector (matrix). */
  DAE_P,
  /** State derivative vector (matrix). Should have same amount of non-zeros as DAEOutput:DAE_RES */
  DAE_YDOT,
  /** Number of arguments. */
  DAE_NUM_IN
};

/// Output arguments of an ODE/DAE residual function
enum DAEOutput{
  /** Right hand side of ODE. Should have same amount of non-zeros as DAEInput:DAE_Y */
  DAE_RES,
  /** Number of arguments. */
  DAE_NUM_OUT
};

/// Input arguments of an ODE/DAE forward integration function (new, not yet ready implementation)
enum DAEFInput{
  /** Differential state */
  DAE_F_XD,
  /** Algebraic state */
  DAE_F_XA,
  /** Parameter vector */
  DAE_F_P,
  /** Explicit time dependence */
  DAE_F_T,
  /** Number of arguments. */
  DAE_F_NUM_IN
};

/// Helper function to create DAE forward integration function input arguments
template<class M>
std::vector<M> daeFIn(const M& xd, const M& xa=M(), const M& p=M(), const M& t=M()){
  M ret[DAE_F_NUM_IN] = {xd,xa,p,t};
  return std::vector<M>(ret,ret+DAE_F_NUM_IN);
}
#ifdef SWIG
%template(daeFIn) daeFIn<SXMatrix>;
%template(daeFIn) daeFIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE forward integration function (new, not yet ready implementation)
enum DAEFOutput{
  /** Right hand side of ODE.*/
  DAE_F_ODE,
  /** Right hand side of algebraic equations.*/
  DAE_F_ALG,
  /** Right hand side of quadratures.*/
  DAE_F_QUAD,
  /** Number of arguments. */
  DAE_F_NUM_OUT
};

/// Helper function to create DAE forward integration function output arguments
template<class M>
std::vector<M> daeFOut(const M& ode, const M& alg=M(), const M& quad=M()){
  M ret[DAE_F_NUM_OUT] = {ode,alg,quad};
  return std::vector<M>(ret,ret+DAE_F_NUM_OUT);
}
#ifdef SWIG
%template(daeFOut) daeFOut<SXMatrix>;
%template(daeFOut) daeFOut<MX>;
#endif //SWIG

/// Input arguments of an ODE/DAE terminal constraint function (new, not yet ready implementation)
enum DAEHInput{
  /** Differential state */
  DAE_H_XD,
  /** Algebraic state */
  DAE_H_XA,
  /** Quadrature state */
  DAE_H_XQ,
  /** Parameter vector */
  DAE_H_P,
  /** Number of arguments. */
  DAE_H_NUM_IN
};

/// Helper function to create ODE/DAE terminal constraint function input arguments
template<class M>
std::vector<M> daeHIn(const M& xd, const M& xa=M(), const M& xq=M(), const M& p=M()){
  M ret[DAE_H_NUM_IN] = {xd,xa,xq,p};
  return std::vector<M>(ret,ret+DAE_H_NUM_IN);
}
#ifdef SWIG
%template(daeHIn) daeHIn<SXMatrix>;
%template(daeHIn) daeHIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE terminal function
enum DAEHOutput{
  /** Initial conditions for the backwards integration, differential states. */
  DAE_H_YD,
  /** Initial conditions for the backwards integration, quadrature states. */
  DAE_H_YA,
  /** Guess for the initial conditions for the backwards integration, algebraic states. */
  DAE_H_YQ,
  /** Number of arguments. */
  DAE_H_NUM_OUT
};

/// Helper function to create DAE terminal constraint function output arguments
template<class M>
std::vector<M> daeHOut(const M& yd, const M& ya=M(), const M& yq=M()){
  M ret[DAE_H_NUM_OUT] = {yd,ya,yq};
  return std::vector<M>(ret,ret+DAE_H_NUM_OUT);
}
#ifdef SWIG
%template(daeHOut) daeHOut<SXMatrix>;
%template(daeHOut) daeHOut<MX>;
#endif //SWIG

/// Input arguments of an ODE/DAE backward integration function (new, not yet ready implementation)
enum DAEGInput{
  /** Forward differential state */
  DAE_G_XD,
  /** Forward algebraic state */
  DAE_G_XA,
  /** Backward differential state */
  DAE_G_YD,
  /** Backward algebraic state */
  DAE_G_YA,
  /** Parameter vector */
  DAE_G_P,
  /** Explicit time dependence */
  DAE_G_T,
  /** Number of arguments. */
  DAE_G_NUM_IN
};

/// Helper function to create ODE/DAE backward integration function input arguments
template<class M>
std::vector<M> daeGIn(const M& xd, const M& xa=M(), const M& yd=M(), const M& ya=M(), const M& p=M(), const M& t=M()){
  M ret[DAE_G_NUM_IN] = {xd,xa,yd,ya,p,t};
  return std::vector<M>(ret,ret+DAE_G_NUM_IN);
}
#ifdef SWIG
%template(daeGIn) daeGIn<SXMatrix>;
%template(daeGIn) daeGIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE backward integration function
enum DAEGOutput{
  /** Right hand side of ODE.*/
  DAE_G_ODE,
  /** Right hand side of algebraic equations.*/
  DAE_G_ALG,
  /** Right hand side of quadratures.*/
  DAE_G_QUAD,
  /** Number of arguments. */
  DAE_G_NUM_OUT
};

/// Helper function to create ODE/DAE backward integration function output arguments
template<class M>
std::vector<M> daeGOut(const M& ode, const M& alg=M(), const M& quad=M()){
  M ret[DAE_G_NUM_OUT] = {ode,alg,quad};
  return std::vector<M>(ret,ret+DAE_G_NUM_OUT);
}
#ifdef SWIG
%template(daeGOut) daeGOut<SXMatrix>;
%template(daeGOut) daeGOut<MX>;
#endif //SWIG

/// Input arguments of an integrator
enum IntegratorInput{
  /** Differential or algebraic state at t0  (dimension nx-by-1) */
  INTEGRATOR_X0, 
  /** Parameters p  (dimension np-by-1) */
  INTEGRATOR_P,  
  /** State derivative at t0  (dimension nx-by-1)
  * Only relevant for implicit intergators.
  * This input may be changed during an IDASIntegrator::evaluate()
  */
  INTEGRATOR_XP0, 
  /** Number of input arguments of an integrator */
  INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum IntegratorOutput{
 /**  State at tf */
 INTEGRATOR_XF, 
 /**  State derivative at tf */
 INTEGRATOR_XPF, 
  /** Number of output arguments of an integrator */
 INTEGRATOR_NUM_OUT
};

/// Input arguments of an integrator (new, not yet ready implementation)
enum NewIntegratorInput{
  /** Differential state at t0 */
  NEW_INTEGRATOR_XD0,
  /** Quadrature state at t0 */
  NEW_INTEGRATOR_XQ0,
  /** Guess for the algebraic state at t0 */
  NEW_INTEGRATOR_XA0,
  /** Parameters p */
  NEW_INTEGRATOR_P,
  /** Number of input arguments of an integrator */
  NEW_INTEGRATOR_NUM_IN};

/// Output arguments of an integrator (new, not yet ready implementation)
enum NewIntegratorOutput{
  /**  Differential state at tf */
  NEW_INTEGRATOR_XDF,
  /**  Quadrature state at tf */
  NEW_INTEGRATOR_XQF,
  /**  Algebraic state at tf */
  NEW_INTEGRATOR_XAF,
  /**  Backward differential state at t0 */
  NEW_INTEGRATOR_YD0,
  /**  Backward quadrature state at t0 */
  NEW_INTEGRATOR_YQ0,
  /**  Backward algebraic state at t0 */
  NEW_INTEGRATOR_YA0,
  /** Number of output arguments of an integrator */
  NEW_INTEGRATOR_NUM_OUT
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
  void reset(int nfdir=0, int nadir=0);

  /// Integrate until a specified time point 
  void integrate(double t_out);

  /// Reset the solver of the adjoint problem and take time to tf
  void resetAdj();

  /// Integrate backwards in time until a specified time point
  void integrateAdj(double t_out);

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Get the DAE
  FX getDAE();
};

} // namespace CasADi

#endif //INTEGRATOR_HPP
