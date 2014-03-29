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
  
  @copydoc DAE_doc
  
*/

namespace CasADi{

// Forward declaration of internal class
class IdasInternal;

/** \brief Interface to IDAS from the Sundials suite.

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
  *   \param g backwards system
  * \copydoc scheme_RDAEInput
  * \copydoc scheme_RDAEOutput
  */
  explicit IdasIntegrator(const Function& f, const Function& g=Function());

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
  static Integrator creator(const Function& f, const Function& g){ return IdasIntegrator(f,g);}
  #ifdef SWIG
  %nocallback;
  #endif
};

} // namespace CasADi

#endif //IDAS_INTEGRATOR_HPP

