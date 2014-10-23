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


#ifndef CASADI_LP_SOLVER_HPP
#define CASADI_LP_SOLVER_HPP

#include "function.hpp"


/** \defgroup LpSolver_doc

    Solves the following linear problem:

    \verbatim
    min          c' x
    x

    subject to
    LBA <= A x <= UBA
    LBX <= x   <= UBX

    with x ( n x 1)
    c ( n x 1 )
    A sparse matrix ( nc x n)
    LBA, UBA dense vector (nc x 1)
    LBX, UBX dense vector (n x 1)

    n: number of decision variables (x)
    nc: number of constraints (A)

    \endverbatim


*/

namespace casadi {

  /// Input arguments of a LP problem [lpIn]
  enum LpSolverInput {
    /// The vector c: dense (n x 1) [c]
    LP_SOLVER_C,
    /// The matrix A: sparse, (nc x n) - product with x must be dense. [a]
    LP_SOLVER_A,
    /// dense, (nc x 1) [lba]
    LP_SOLVER_LBA,
    /// dense, (nc x 1) [uba]
    LP_SOLVER_UBA,
    /// dense, (n x 1) [lbx]
    LP_SOLVER_LBX,
    /// dense, (n x 1) [ubx]
    LP_SOLVER_UBX,
    LP_SOLVER_NUM_IN};

  /// Output arguments of an LP Solver [lpOut]
  enum LpSolverOutput {
    /// The primal solution [x]
    LP_SOLVER_X,
    /// The optimal cost [cost]
    LP_SOLVER_COST,
    /// The dual solution corresponding to linear bounds [lam_a]
    LP_SOLVER_LAM_A,
    /// The dual solution corresponding to simple bounds [lam_x]
    LP_SOLVER_LAM_X,
    LP_SOLVER_NUM_OUT};

  /// Structure specification of an LP [lpStruct]
  enum LPStruct {
    /// The matrix A: sparse [a]
    LP_STRUCT_A,
    LP_STRUCT_NUM};

  // Forward declaration of internal class
  class LpSolverInternal;

  /** \brief LpSolver


      @copydoc LpSolver_doc

      \generalsection{LpSolver}
      \pluginssection{LpSolver}

      \author Joris Gillis
      \date 2013
  */
  class CASADI_CORE_EXPORT LpSolver : public Function {
  public:

    /// Default constructor
    LpSolver();

    /** \brief Constructor
     *  \param name \pluginargument{LpSolver}
     *  \param st \structargument{LP}
     */
    LpSolver(const std::string& name, const LPStructure& st);

    /// Access functions of the node
    LpSolverInternal* operator->();
    const LpSolverInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);
  };

} // namespace casadi

#endif // CASADI_LP_SOLVER_HPP

