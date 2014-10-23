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


#ifndef CASADI_DLE_SOLVER_HPP
#define CASADI_DLE_SOLVER_HPP

#include "function.hpp"

/** \defgroup DLE_doc Discrete periodic Lyapunov Equation solver


    Given matrices \f$A\f$ and symmetric \f$V\f$

    \verbatim
    A in R^(n x n)
    V in S^n
    \endverbatim

    finds \f$P\f$ that satisfies:

    \verbatim
    P = A P A' + V
    \endverbatim


*/
namespace casadi {

  /// Input arguments of a \e dle solver [dleIn]
  enum DLEInput {
    /// A matrix [a]
    DLE_A,
    /// V matrix [v]
    DLE_V,
    DLE_NUM_IN
  };

  /// Output arguments of a \e dle solver [dleOut]
  enum DLEOutput {
    /// P matrix [p]
    DLE_P,
    DLE_NUM_OUT
  };

  /// Structure specification of a DLE [dleStruct]
  enum DleStruct {
    /// The matrix A [a]
    Dle_STRUCT_A,
    /// The matrix V [v]
    Dle_STRUCT_V,
    Dle_STRUCT_NUM};

  /// Forward declaration of internal class
  class DleInternal;

  /**  \brief Base class for Discrete Lyapunov Equation Solvers

     @copydoc DLE_doc

      \generalsection{DleSolver}
      \pluginssection{DleSolver}

       \author Joris Gillis
      \date 2014

  */
  class CASADI_CORE_EXPORT DleSolver : public Function {
  public:
    /// Default constructor
    DleSolver();

    /// Clone
    DleSolver clone() const;

    /** \brief DleSolver solver factory
    * \param name \pluginargument{DleSolver}
    * \param st \structargument{Dle}
    */
    DleSolver(const std::string& name, const DleStructure& st);

    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;

    /// Access functions of the node
    DleInternal* operator->();

    /// Access functions of the node
    const DleInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Get the resulting sparsity
    static Sparsity getSparsity(const DleStructure& st);
  };

} // namespace casadi

#endif // CASADI_DLE_SOLVER_HPP
