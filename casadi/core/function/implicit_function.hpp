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


#ifndef CASADI_IMPLICIT_FUNCTION_HPP
#define CASADI_IMPLICIT_FUNCTION_HPP

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

      \generalsection{ImplicitFunction}
      \pluginssection{ImplicitFunction}

     \author Joel Andersson
     \date 2011
  */
  class CASADI_CORE_EXPORT ImplicitFunction : public Function {
  public:

    /** \brief  Default constructor */
    ImplicitFunction();

    /** \brief  Create an implicit function solver
     * \param name \pluginargument{ImplicitFunction}
     * \param f Function mapping from (n+1) inputs to 1 output
     *
     */
    ImplicitFunction(const std::string& name, const Function& f,
                     const Function& jac=Function(),
                     const LinearSolver& linsol=LinearSolver());

    /// Access functions of the node
    ImplicitFunctionInternal* operator->();

    /// Const access functions of the node
    const ImplicitFunctionInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Access F
    Function& getF();

    /// Access Jacobian
    Function& getJac();

    /// Access linear solver
    LinearSolver& getLinsol();
  };

} // namespace casadi

#endif // CASADI_IMPLICIT_FUNCTION_HPP

