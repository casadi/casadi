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


#ifndef CASADI_CLE_SOLVER_HPP
#define CASADI_CLE_SOLVER_HPP

#include "function.hpp"

/** \defgroup CLE_doc Discrete periodic Lyapunov Equation solver


    Given matrices \f$A\f$ and symmetric \f$V\f$

    \verbatim
    A in R^(n x n)
    V in S^n
    \endverbatim

    finds \f$P\f$ that satisfies:

    \verbatim
    0 = A P  + P A' + V
    \endverbatim


*/
namespace casadi {

  /// Input arguments of a \e cle solver [cleIn]
  enum CLEInput {
    /// A matrix [a]
    CLE_A,
    /// V matrix [v]
    CLE_V,
    CLE_NUM_IN
  };

  /// Output arguments of a \e cle solver [cleOut]
  enum CLEOutput {
    /// Lyapunov matrix [p]
    CLE_P,
    /// Number of arguments.
    CLE_NUM_OUT
  };

  /// Structure specification of a CLE [cleStruct]
  enum CleStruct {
    /// The matrix A [a]
    Cle_STRUCT_A,
    /// The matrix V [v]
    Cle_STRUCT_V,
    /// The matrix C (defaults to unity) [c]
    Cle_STRUCT_C,
    Cle_STRUCT_NUM};


  /// Forward declaration of internal class
  class CleInternal;

  /**  \brief Base class for Discrete Lyapunov Equation Solvers

     @copydoc CLE_doc

      \generalsection{CleSolver}
      \pluginssection{CleSolver}

       \author Joris Gillis
      \date 2014

  */
  class CASADI_CORE_EXPORT CleSolver : public Function {
  public:
    /// Default constructor
    CleSolver();

    /// Clone
    CleSolver clone() const;

    /** \brief CleSolver solver factory
    * \param name \pluginargument{CleSolver}
    * \param st \structargument{Cle}
    */
    CleSolver(const std::string& name,
               const CleStructure& st);

    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;

    /// Access functions of the node
    CleInternal* operator->();

    /// Access functions of the node
    const CleInternal* operator->() const;

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

#endif // CASADI_CLE_SOLVER_HPP
