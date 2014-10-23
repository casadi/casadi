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


#ifndef CASADI_LR_DLE_SOLVER_HPP
#define CASADI_LR_DLE_SOLVER_HPP

#include "function.hpp"

/** \defgroup LR_DLE_doc Low Rank Discrete periodic Lyapunov Equation solver


    \verbatim
    A in R^(n x n)
    V in S^m
    C in R^(n x m)
    Hi in R^(n x Hsi)
    \endverbatim

    finds \f$P\f$ that satisfies:

    \verbatim
    P = A P A' + C V C'
    \endverbatim

    and outputs

    Yi = Hi^T P Hi


*/
namespace casadi {

  /// Input arguments of a \e dle solver [lrdleIn]
  enum LR_DLEInput {
    /// A matrix [a]
    LR_DLE_A,
    /// V matrix [v]
    LR_DLE_V,
    /// C matrix [c]
    LR_DLE_C,
    /// H matrix: horizontal stack of all Hi [h]
    LR_DLE_H,
    LR_DLE_NUM_IN
  };

  /// Output arguments of a \e dle solver [lrdleOut]
  enum LR_DLEOutput {
    /// Y matrix, blkdiag form [y]
    LR_DLE_Y,
    LR_DLE_NUM_OUT
  };

  /// Structure specification of a DLE [lrdleStruct]
  enum LrDleStruct {
    /// The matrix A [a]
    LR_DLE_STRUCT_A,
    /// The matrix V [v]
    LR_DLE_STRUCT_V,
    /// The matrix C (defaults to unity) [c]
    LR_DLE_STRUCT_C,
    /// H matrix: horizontal stack of all Hi [h]
    LR_DLE_STRUCT_H,
    LR_DLE_STRUCT_NUM};

  /// Forward declaration of internal class
  class LrDleInternal;

  /**  \brief Base class for Low-rank Discrete Lyapunov Equation Solvers

     @copydoc LR_DLE_doc

      \generalsection{LrDleSolver}
      \pluginssection{LrDleSolver}

       \author Joris Gillis
      \date 2014

  */
  class CASADI_CORE_EXPORT LrDleSolver : public Function {
  public:
    /// Default constructor
    LrDleSolver();

    /// Clone
    LrDleSolver clone() const;

    /** \brief LrDleSolver solver factory
    * \param name \pluginargument{LrDleSolver}
    * \param st \structargument{LrDle}
    * \param Hs Column-sizes of H_i
    */
    LrDleSolver(const std::string& name, const LrDleStructure& st,
      const std::vector<int> &Hs=std::vector<int>());

    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;

    /// Access functions of the node
    LrDleInternal* operator->();

    /// Access functions of the node
    const LrDleInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Get the resulting sparsity
    static Sparsity getSparsity(const LrDleStructure& st,
      const std::vector<int> &Hs=std::vector<int>());
  };

} // namespace casadi

#endif // CASADI_LR_DLE_SOLVER_HPP
