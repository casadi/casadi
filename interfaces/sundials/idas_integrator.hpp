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

#ifndef IDAS_INTEGRATOR_HPP
#define IDAS_INTEGRATOR_HPP

#include "sundials_integrator.hpp"

/**
* \defgroup IdasIntegrator_doc
  Solves an initial value problem in differential-algebraic equations of the form:
  
  
  Creates an integrator instance which solves initial value problems in differential-algebraic equations
  of the form:
  
  \verbatim
  f(t,y,der(y),z,p) == 0
  der(q) = g(t,y,z,p)
  \endverbatim
  
  The DAE thus consists of a fully implicit part (f) and an explicit quadrature part (g). In the same way,
  the state vector is also composed of two parts, the differential states and the quadrature states,
  i.e. x = [y,q]
*/

namespace CasADi{

// Forward declaration of internal class
class IdasInternal;

/// Input arguments of a jacobian function: J = [df/dx + cj*df/dxdot, df/dz]
enum JACInput{JAC_T, JAC_X, JAC_Z, JAC_XDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/// Output arguments of an DAE residual function
enum JACOutput{JAC_J, JAC_NUM_OUT};

/** Interface to IDAS from the Sundials suite.

   @copydoc IdasIntegrator_doc
  
   \author Joel Andersson
   \date 2010
*/

class IdasIntegrator : public SundialsIntegrator{
public:

  /// Default constructor
  IdasIntegrator();
  
  /// Create an integrator for a fully implicit DAE with quadrature states (nz is the number of states not to be included in the state vector)
  
    
  /** \brief  Create an integrator for a fully implicit DAE with quadrature states 
  * (nz is the number of states not to be included in the state vector)
  *   \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  *
  */
  explicit IdasIntegrator(const FX& f, const FX& g=FX());

  /// Access functions of the node
  IdasInternal* operator->();

  /// Const access functions of the node
  const IdasInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Correct the initial value for yp and z after resetting the solver
  void correctInitialConditions();

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static Integrator creator(const FX& f, const FX& g){ return IdasIntegrator(f,g);}
  #ifdef SWIG
  %nocallback;
  #endif
};

} // namespace CasADi

#endif //IDAS_INTEGRATOR_HPP

