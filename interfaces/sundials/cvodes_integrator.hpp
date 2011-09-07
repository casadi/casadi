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

#ifndef CVODES_INTEGRATOR_HPP
#define CVODES_INTEGRATOR_HPP

#include "casadi/fx/integrator.hpp"




namespace CasADi{
namespace Sundials{
  
/// Input arguments of the Jacobian in the nonlinear iteration: M = 1 - gamma*df/dy
enum MInput{M_T, M_Y, M_P, M_GAMMA, M_NUM_IN};

/// Output arguments of the Jacobian function
enum MOutput{M_M, M_NUM_OUT};
  
// Forward declaration of internal class 
class CVodesInternal;

/** Function that integrates the ODE:

  @copydoc ODE_doc
  
  A call to evaluate will integrate to the end.
  
  You can retrieve the entire state trajectory as follows, after the evaluate call: 
  Call reset. Then call integrate(t_i) and getOuput for a series of times t_i.
  

*/
class CVodesIntegrator : public Integrator{
public:

  /** \brief  Default constructor */
  CVodesIntegrator();
  
  /** \brief  Create an integrator for explicit ODEs
  *   \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  *
  */
  explicit CVodesIntegrator(const FX& f, const FX& q=FX());
  
  /** \brief  Access functions of the node */
  CVodesInternal* operator->();
  const CVodesInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static Integrator creator(const FX& f, const FX& q){ return CVodesIntegrator(f,q);}
  #ifdef SWIG
  %nocallback;
  #endif
};


} // namespace Sundials
} // namespace CasADi

#endif //CVODES_INTEGRATOR_HPP

