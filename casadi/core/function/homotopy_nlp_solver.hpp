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


#ifndef CASADI_HOMOTOPY_NLP_SOLVER_HPP
#define CASADI_HOMOTOPY_NLP_SOLVER_HPP

#include "function.hpp"


/** \defgroup HomotopyNlpSolver_doc

  Solves the following parametric nonlinear program (NLP):
  \verbatim
  min          F(x, p, tau)
   x

  subject to
              LBX <=   x    <= UBX
              LBG <= G(x, p) <= UBG
                         p  == P

      nx: number of decision variables
      ng: number of constraints
      np: number of parameters
  \endverbatim

  In a homotopy from tau = 0 to tau = 1.

*/

namespace casadi {

  /// Input arguments of an Homotopy NLP function [hnlpIn]
  enum HNLPInput {
    /// Decision variable [x]
    HNL_X,
    /// Fixed parameter [p]
    HNL_P,
    /// Homotopy parameter [tau]
    HNL_TAU,
    /// Number of NLP inputs
    HNL_NUM_IN
  };

  class HomotopyNLPInternal;

  /** \brief Base class for Homotopy NLP Solvers

      @copydoc HomotopyNlpSolver_doc

      \generalsection{HomotopyNlpSolver}
      \pluginssection{HomotopyNlpSolver}

      \author Joris Gillis
      \date 2014
  */
  class CASADI_CORE_EXPORT HomotopyNlpSolver : public Function {
  public:

    /// Default constructor
    HomotopyNlpSolver();

    /// Access functions of the node
    HomotopyNLPInternal* operator->();
    const HomotopyNLPInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);
  };

} // namespace casadi

#endif // CASADI_HOMOTOPY_NLP_SOLVER_HPP

