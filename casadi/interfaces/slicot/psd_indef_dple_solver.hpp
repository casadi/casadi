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


#ifndef CASADI_PSD_INDEF_DPLE_SOLVER_HPP
#define CASADI_PSD_INDEF_DPLE_SOLVER_HPP

#include "../../core/function/dple_solver.hpp"
#include <casadi/interfaces/slicot/casadi_dplesolver_slicot_export.h>

namespace casadi {

  /// Forward declaration of internal class
  class PsdIndefDpleInternal;

  /** \brief An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT
   *
   * @copydoc DPLE_doc

       Uses Periodic Schur Decomposition ('psd') and does not assume positive definiteness.
       Based on Periodic Lyapunov equations: some applications and new algorithms.
       Int. J. Control, vol. 67, pp. 69-87, 1997.

       \author Joris Gillis
      \date 2014

  */
  class CASADI_DPLESOLVER_SLICOT_EXPORT PsdIndefDpleSolver : public DpleSolver {
  public:
    /// Default constructor
    PsdIndefDpleSolver();

    /** \brief  Constructor
     *  \param[in] A  List of sparsities of A_i
     *  \param[in] V  List of sparsities of V_i
     */
    explicit PsdIndefDpleSolver(const std::vector< Sparsity > & A,
                                const std::vector< Sparsity > &V);

    /// Access functions of the node
    PsdIndefDpleInternal* operator->();

    /// Access functions of the node
    const PsdIndefDpleInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Static creator function
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static DpleSolver creator(const std::vector< Sparsity > & A, const std::vector< Sparsity > &V)
    { return PsdIndefDpleSolver(A, V); }
    #ifdef SWIG
    %nocallback;
    #endif

  };

} // namespace casadi

#endif // CASADI_DPLE_SOLVER_HPP
