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

// To be able to maintain both Python and Doxygen documentation
#ifdef SWIG
#define SWIGDOC(x) \
%feature("autodoc", x " Consult the C++ doxygen information for details.");
#else
#define SWIGDOC(x)
#endif

namespace CasADi{
namespace Sundials{

// Forward declaration of internal class
class IdasInternal;

/// Input arguments of an DAE residual function
enum DAEInput{DAE_T, DAE_Y, DAE_YDOT, DAE_P, DAE_NUM_IN};

// Input arguments of an DAE residual function (new version)
//enum DAEInput{DAE_T, DAE_Y, DAE_YDOT, DAE_Z, DAE_P, DAE_NUM_IN};

/// Output arguments of an DAE residual function
enum DAEOutput{DAE_RES, DAE_NUM_OUT};

/// Input arguments of a jacobian function: J = df/dy + cj*df/dydot
enum JACInput{JAC_T, JAC_Y, JAC_YDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/// Output arguments of an DAE residual function
enum JACOutput{JAC_J, JAC_NUM_OUT};

/** Interface to IDAS from the Sundials suite.

  Creates an integrator instance which solves initial value problems in differential-algebraic equations
  of the form:
  
  f(t,y,der(y),z,p) == 0
  der(q) = g(t,y,z,p)

  The DAE thus consists of a fully implicit part (f) and an explicit quadrature part (g). In the same way,
  the state vector is also composed of two parts, the differential states and the quadrature states,
  i.e. x = [y,q]
  
   \author Joel Andersson
   \date 2010
*/

SWIGDOC("Interface to IDAS from the Sundials suite.");
class IdasIntegrator : public Integrator{
public:

  /// Default constructor
  IdasIntegrator();
  
  /// Create an integrator for a fully implicit DAE with quadrature states
  explicit IdasIntegrator(const FX& f, const FX& q=FX());

  /// Access functions of the node
  IdasInternal* operator->();

  /// Const access functions of the node
  const IdasInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Generate a new integrator integrating the forward sensitivity augmented ODE/DAE
  IdasIntegrator jac(int iind=0, int oind=0);
};


} // namespace Sundials
} // namespace CasADi

#endif //IDAS_INTEGRATOR_HPP

