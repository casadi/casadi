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

      Mathematically, the equation:

      F(z, x1, x2, ..., xn) == 0

      where d_F/dz is invertible, implicitly defines the equation:

      z := G(x1, x2, ..., xn)

      In CasADi, F is a Function.
      The first input presents the variables that need to be solved for.
      The first output is the residual that needs to attain zero.
      Optional remaining outputs can be supplied; they are output expressions.
      
      In pseudo-code, we can write:
      
      G* = ImplicitFunction('solver',F)
      
      Here, G* is a Function with one extra input over the pure mathematical G:
      
      z := G*(z0, x1, x2, ..., xn)

      The first input to the ImplicitFunction is the intial guess for z.

  */
  /** \brief  Abstract base class for the implicit function classes

     @copydoc ImplicitFunction_doc

      \generalsection{ImplicitFunction}
      \pluginssection{ImplicitFunction}

     \author Joel Andersson
     \date 2011
  */
  class CASADI_EXPORT ImplicitFunction : public Function {
  public:

    /** \brief  Default constructor */
    ImplicitFunction();

    /** \brief  Create an implicit function solver (new syntax, includes initialization)
     * \param solver \pluginargument{ImplicitFunction}
     * \param f Function where one of the inputs (by default the first) is an unknown and
     *        one of the outputs (by default the first) is a residual.
     */
    ImplicitFunction(const std::string& name, const std::string& solver,
                     const Function& f, const Dict& opts=Dict());

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Create an implicit function solver, no initialization
     * \param solver \pluginargument{ImplicitFunction}
     * \param f Function where one of the inputs (by default the first) is an unknown and
     *        one of the outputs (by default the first) is a residual.
     */
    ImplicitFunction(const std::string& solver, const Function& f,
                     const Function& jac=Function(),
                     const LinearSolver& linsol=LinearSolver());
#endif // WITH_DEPRECATED_FEATURES

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

