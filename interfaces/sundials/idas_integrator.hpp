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

#include "casadi/fx/integrator.hpp"

/** Integrator where the differential equation is assumed to be in the form:

  index-1 DAE:
  f(t,y,der(y),p) == 0

  pure quadratures:
  der(q) = g(t,y,p)
  
  The state vector is x = [x,q]
  
*/

namespace CasADi{
namespace Sundials{

/** \brief  Forward declaration of internal class */
class IdasInternal;

/** \brief  Input arguments of a jacobian function: J = df/dy + cj*df/dydot */
enum JACInput{JAC_T, JAC_Y, JAC_YDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/** \brief  Output arguments of an DAE residual function */
enum JACOutput{JAC_J, JAC_NUM_OUT};

/** \brief  Public class */
class IdasIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  IdasIntegrator();
  
  /** \brief  Create an integrator for explicit ODEs */
  explicit IdasIntegrator(const FX& f, const FX& q=FX());

  /** \brief  Access functions of the node */
  IdasInternal* operator->();
  const IdasInternal* operator->() const;
  
  /** \brief  Set linear solver */
  void setLinearSolver(const FX& jac, const FX& linsol);
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

};


} // namespace Sundials
} // namespace CasADi

#endif //IDAS_INTEGRATOR_HPP

