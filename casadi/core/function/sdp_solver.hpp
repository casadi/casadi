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


#ifndef CASADI_SDP_SOLVER_HPP
#define CASADI_SDP_SOLVER_HPP

#include "function.hpp"


/** \defgroup SdpSolver_doc

    Solves an SDP problem in standard form.
    See http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf

    Primal:

    \verbatim
    min          c' x
    x
    subject to
    P = Sum_i^m F_i x_i - G
    P negative semidefinite

    LBA <= A x <= UBA
    LBX <= x   <= UBX

    with x ( n x 1)
    c ( n x 1 )
    G, F_i  sparse symmetric (m x m)
    X dense symmetric ( m x m )
    A sparse matrix ( nc x n)
    LBA, UBA dense vector (nc x 1)
    LBX, UBX dense vector (n x 1)

    \endverbatim

    This formulation is chosen as primal, because it does not call for a large decision variable space.

    Dual:

    \verbatim
    max          trace(G Y)
    Y

    subject to
    trace(F_i Y) = c_i
    Y positive semidefinite

    with Y dense symmetric ( m x m)

    \endverbatim

    On generality: you might have formulation with block partitioning:

    Primal:

    \verbatim
    min          c' x
    x
    subject to
    Pj = Sum_i^m F_ij x_i - gj   for all j
    Pj negative semidefinite   for all j

    with x ( n x 1)
    c ( n x 1 )
    G, F_i  sparse symmetric (m x m)
    X dense symmetric ( m x m )

    \endverbatim

    Dual:
    \verbatim
    max          Sum_j trace(Gj Yj)
    Yj

    subject to
    Sum_j trace(F_ij Yj) = c_i   for all j
    Yj positive semidefinite     for all j

    with Y dense symmetric ( m x m)

    \endverbatim

    You can cast this into the standard form with:
    G  = blkdiag(Gj for all j)
    Fi = blkdiag(F_ij for all j)

    Implementations of SdpSolver are encouraged to exploit this block structure.

*/

namespace casadi {

  /// Input arguments of a SDP problem [sdpIn]
  enum SDPInput {
    /// The horizontal stack of all matrices F_i: ( m x nm) [f]
    SDP_SOLVER_F,
    /// The vector c: ( n x 1) [c]
    SDP_SOLVER_C,
    /// The matrix G: ( m x m) [g]
    SDP_SOLVER_G,
    /// The matrix A: ( nc x n) [a]
    SDP_SOLVER_A,
    /// Lower bounds on Ax ( nc x 1) [lba]
    SDP_SOLVER_LBA,
    /// Upper bounds on Ax  ( nc x 1) [uba]
    SDP_SOLVER_UBA,
    /// Lower bounds on x ( n x 1 ) [lbx]
    SDP_SOLVER_LBX,
    /// Upper bounds on x ( n x 1 ) [ubx]
    SDP_SOLVER_UBX,
    SDP_SOLVER_NUM_IN};

  /// Output arguments of an SDP Solver [sdpOut]
  enum SDPOutput {
    /// The primal solution (n x 1) - may be used as initial guess [x]
    SDP_SOLVER_X,
    /// The solution P (m x m) - may be used as initial guess [p]
    SDP_SOLVER_P,
    /// The dual solution (m x m) - may be used as initial guess [dual]
    SDP_SOLVER_DUAL,
    /// The primal optimal cost (1 x 1) [cost]
    SDP_SOLVER_COST,
    /// The dual optimal cost (1 x 1) [dual_cost]
    SDP_SOLVER_DUAL_COST,
    /// The dual solution corresponding to the linear constraints  (nc x 1) [lam_a]
    SDP_SOLVER_LAM_A,
    /// The dual solution corresponding to simple bounds  (n x 1) [lam_x]
    SDP_SOLVER_LAM_X,
    SDP_SOLVER_NUM_OUT};

  /// Structure specification of an SDP [sdpStruct]
  enum SDPStruct {
    /// The horizontal stack of all matrices F_i: ( m x nm) [f]
    SDP_STRUCT_F,
    /// The matrix G: ( m x m) [g]
    SDP_STRUCT_G,
    /// The matrix A: ( nc x n) [a]
    SDP_STRUCT_A,
    SDP_STRUCT_NUM};

  // Forward declaration of internal class
  class SdpSolverInternal;

  /** \brief SdpSolver


      @copydoc SdpSolver_doc

      \generalsection{SdpSolver}
      \pluginssection{SdpSolver}

      \author Joel Andersson
      \date 2010
  */
  class CASADI_CORE_EXPORT SdpSolver : public Function {
  public:

    /// Default constructor
    SdpSolver();

    /** \brief Constructor
     *  \param name \pluginargument{SdpSolver}
     *  \param st \structargument{SDP}
     */
    SdpSolver(const std::string& name, const SDPStructure& st);

    /// Access functions of the node
    SdpSolverInternal* operator->();
    const SdpSolverInternal* operator->() const;

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Set options that make the SDP solver more suitable for solving SOCPs
    void setSOCPOptions();

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_SDP_SOLVER_HPP
