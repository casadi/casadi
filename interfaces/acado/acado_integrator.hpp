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

#ifndef ACADO_INTEGRATOR_HPP
#define ACADO_INTEGRATOR_HPP

#include "casadi/fx/integrator.hpp"

/**
* \defgroup AcadoIntegrator_doc
  Creates an integrator instance which solves initial value problems in differential-algebraic equations
  of the semi-explicit form:
  
  \verbatim
  der(xd) = fd(t,[xd,xa],p)
      0   = fa(t,[xd,xa],p)
  \endverbatim
  
  The DAE thus consists of an differential part (fd) and an algebraic part (fa). In the same way,
  the state vector is also composed of two parts, the differential states xd and the algebraic states
  xa. The complete state is thus x := [xd,xa]
*/

namespace CasADi{

// Forward declaration of internal class
class AcadoIntegratorInternal;

/** Interface to ACADO Integrator from the ACADO Toolkit.

   @copydoc AcadoIntegrator_doc
  
   \author Joel Andersson
   \date 2011
*/

class AcadoIntegrator : public Integrator{
public:

  /// Default constructor
  AcadoIntegrator();
  
  /// Create an integrator for a fully implicit DAE with quadrature states (nz is the number of states not to be included in the state vector)
  explicit AcadoIntegrator(const FX& f, const FX& q=FX());

  /// Access functions of the node
  AcadoIntegratorInternal* operator->();

  /// Const access functions of the node
  const AcadoIntegratorInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static Integrator creator(const FX& f, const FX& q){ return AcadoIntegrator(f,q); }
  #ifdef SWIG
  %nocallback;
  #endif
};


} // namespace CasADi

#endif //ACADO_INTEGRATOR_HPP

