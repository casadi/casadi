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
  
  NOTE: The ODE/DAE initial-value problem formulation in CasADi is being replaced with a new 
  two-point boundary value problem with the differential equation given as a semi-explicit 
  DAE with quadrature states:
  \verbatim
  Initial conditions at t=t0
    x(t0)  = x_0
    q(t0)  = 0
  
  Forward integration from t=t0 to t=tf
    der(x) = fx(x,z,p,t)           Forward ODE
         0 = fz(x,z,p,t)           Forward algebraic equations
    der(q) = fq(x,z,p,t)           Forward quadratures
  
  Terminal conditions at t=tf
    rx(tf)  = h(x(tf),q(tf),p)
    rq(tf)  = 0
  
  Backward integration from t=tf to t=t0
    der(rx) = gx(x,z,rx,rz,p,t)   Backward ODE
          0 = gz(x,z,rx,rz,p,t)   Backward algebraic equations
    der(rq) = gq(x,z,rx,rz,p,t)   Backward quadratures

  Where we assume that both the forward and backwards integrations are index-1
  (i.e. fz/dz and dgz/drz are invertible) and furthermore that 
  gx, gz and gq have a linear dependency on rx and rz.
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

/// Input arguments of an ODE/DAE forward integration function
enum NEW_DAEInput{
  /** Differential state */
  NEW_DAE_X,
  /** Algebraic state */
  NEW_DAE_Z,
  /** Parameter vector */
  NEW_DAE_P,
  /** Explicit time dependence */
  NEW_DAE_T,
  /** Number of arguments. */
  NEW_DAE_NUM_IN
};

/// Helper function to create ODE/DAE forward integration function input arguments
template<class M>
std::vector<M> daeIn(const M& x, const M& z=M(), const M& p=M(), const M& t=M()){
  M ret[NEW_DAE_NUM_IN] = {x,z,p,t};
  return std::vector<M>(ret,ret+NEW_DAE_NUM_IN);
}
#ifdef SWIG
%template(daeIn) daeIn<SXMatrix>;
%template(daeIn) daeIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE forward integration function
enum NEW_DAEOutput{
  /** Right hand side of ODE.*/
  NEW_DAE_ODE,
  /** Right hand side of algebraic equations.*/
  NEW_DAE_ALG,
  /** Right hand side of quadratures.*/
  NEW_DAE_QUAD,
  /** Number of arguments. */
  NEW_DAE_NUM_OUT
};

/// Helper function to create ODE/DAE forward integration function output arguments
template<class M>
std::vector<M> daeOut(const M& ode, const M& alg=M(), const M& quad=M()){
  M ret[NEW_DAE_NUM_OUT] = {ode,alg,quad};
  return std::vector<M>(ret,ret+NEW_DAE_NUM_OUT);
}
#ifdef SWIG
%template(daeOut) daeOut<SXMatrix>;
%template(daeOut) daeOut<MX>;
#endif //SWIG

/// Input arguments of an ODE/DAE terminal constraint function
enum TermInput{
  /** Differential state */
  TERM_X,
  /** Quadrature state */
  TERM_Q,
  /** Parameter vector */
  TERM_P,
  /** Number of arguments. */
  TERM_NUM_IN
};

/// Helper function to create ODE/DAE terminal constraint function input arguments
template<class M>
std::vector<M> termIn(const M& x, const M& z=M(), const M& q=M(), const M& p=M()){
  M ret[TERM_NUM_IN] = {x,z,q,p};
  return std::vector<M>(ret,ret+TERM_NUM_IN);
}
#ifdef SWIG
%template(termIn) termIn<SXMatrix>;
%template(termIn) termIn<MX>;
#endif //SWIG

/// Output arguments of an ODE/DAE terminal function
enum TermOutput{
  /** Initial conditions for the backwards integration, differential states. */
  TERM_RX,
  /** Number of arguments. */
  TERM_NUM_OUT
};

/// Helper function to create ODE/DAE terminal constraint function output arguments
template<class M>
std::vector<M> termOut(const M& rx, const M& rz=M(), const M& rq=M()){
  M ret[TERM_NUM_OUT] = {rx,rz,rq};
  return std::vector<M>(ret,ret+TERM_NUM_OUT);
}
#ifdef SWIG
%template(termOut) termOut<SXMatrix>;
%template(termOut) termOut<MX>;
#endif //SWIG

/// Input arguments of an ODE/DAE backward integration function 
enum RDAEInput{
  /** Forward differential state */
  RDAE_X,
  /** Forward algebraic state */
  RDAE_Z,
  /** Backward differential state */
  RDAE_RX,
  /** Backward algebraic state */
  RDAE_RZ,
  /** Parameter vector */
  RDAE_P,
  /** Explicit time dependence */
  RDAE_T,
  /** Number of arguments. */
  RDAE_NUM_IN
};

/// Helper function to create ODE/DAE backward integration function input arguments
template<class M>
std::vector<M> rdaeIn(const M& x, const M& z=M(), const M& rx=M(), const M& rz=M(), const M& p=M(), const M& t=M()){
  M ret[RDAE_NUM_IN] = {x,z,rx,rz,p,t};
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

/// Input arguments of an integrator 
enum NewIntegratorInput{
  /** Differential state at t0 */
  NEW_INTEGRATOR_X0,
  /** Parameters p */
  NEW_INTEGRATOR_P,
  /** Number of input arguments of an integrator */
  NEW_INTEGRATOR_NUM_IN};

/// Output arguments of an integrator
enum NewIntegratorOutput{
  /**  Differential state at tf */
  NEW_INTEGRATOR_XF,
  /**  Quadrature state at tf */
  NEW_INTEGRATOR_QF,
  /**  Backward differential state at t0 */
  NEW_INTEGRATOR_RX0,
  /**  Backward quadrature state at t0 */
  NEW_INTEGRATOR_RQ0,
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
