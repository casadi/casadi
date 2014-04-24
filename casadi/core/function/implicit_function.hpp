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

#ifndef IMPLICIT_FUNCTION_HPP
#define IMPLICIT_FUNCTION_HPP

#include "function.hpp"
#include "linear_solver.hpp"

namespace casadi {
  // Forward declaration of internal class
  class ImplicitFunctionInternal;

  /**
      \defgroup ImplicitFunction_doc

      The equation:

      F(z, x1, x2, ..., xn) == 0

      where d_F/dz is invertible, implicitly defines the equation:

      z := G(x1, x2, ..., xn)



      F should be an Function mapping from (n+1) inputs to m outputs.
      The first output is the residual that should be zero.

      ImplicitFunction (G) is an Function mapping from n inputs to m outputs.
      n may be zero.
      The first output is the solved for z.

      You can provide an initial guess for z by setting output(0) of ImplicitFunction.


  */
  /** \brief  Abstract base class for the implicit function classes

     @copydoc ImplicitFunction_doc

     \author Joel Andersson
     \date 2011
  */
  class CASADI_CORE_EXPORT ImplicitFunction : public Function {
  public:

    /// Access functions of the node
    ImplicitFunctionInternal* operator->();

    /// Const access functions of the node
    const ImplicitFunctionInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Access F
    Function& getF();

    /// Access Jacobian
    Function& getJac();

    /// Access linear solver
    LinearSolver& getLinsol();
  };

} // namespace casadi

#endif //IMPLICIT_FUNCTION_HPP

