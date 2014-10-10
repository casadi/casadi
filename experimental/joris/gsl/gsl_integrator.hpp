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


#ifndef GSL_INTEGRATOR_HPP
#define GSL_INTEGRATOR_HPP

#include "core/function/integrator.hpp"


// http://www.network-theory.co.uk/docs/gslref/OrdinaryDifferentialEquations.html

namespace casadi {
namespace GSL {
  
  
// Forward declaration of internal class 
class GslInternal;

/** 
 @copydoc ODE_doc
  
  A call to evaluate will integrate to the end.
  
  You can retrieve the entire state trajectory as follows, after the evaluate call: 
  Call reset. Then call integrate(t_i) and getOuput for a series of times t_i.
  

*/
class GslIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  GslIntegrator();
  
  /** \brief  Create an integrator for explicit ODEs
  *   \param f casadi::Function mapping from casadi::ODEInput to casadi::ODEOutput.
  */
  explicit GslIntegrator(const Function& f, const Function& q=Function());
  
  /** \brief  Access functions of the node */
  GslInternal* operator->();
  const GslInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
};


} // namespace GSL
} // namespace casadi

#endif //GSL_INTEGRATOR_HPP

