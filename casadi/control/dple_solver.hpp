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

#ifndef CASADI_DPLE_SOLVER_HPP
#define CASADI_DPLE_SOLVER_HPP

#include "../core/function/function.hpp"

#include <casadi/control/casadi_control_export.h>

/** \defgroup DPLE_doc Discrete periodic Lyapunov Equation solver


    Given matrices \f$A_k\f$ and symmetric \f$V_k,  k = 0..K-1\f$

    \verbatim
    A_k in R^(n x n)
    V_k in R^n
    \endverbatim

    provides all of \f$P_k\f$ that satisfy:

    \verbatim
    P_0 = A_(K-1)*P_(K-1)*A_(K-1)' + V_k
    P_k+1 = A_k*P_k*A_k' + V_k  for k = 1..K-1
    \endverbatim


*/
namespace casadi {

  /// Input arguments of a \e dple solver [dpleIn]
  enum DPLEInput {
    /// A matrices (horzcat when const_dim, blkdiag otherwise) [a]
    DPLE_A,
    /// V matrices (horzcat when const_dim, blkdiag otherwise) [v]
    DPLE_V,
    DPLE_NUM_IN
  };

  /// Output arguments of a \e dple solver [dpleOut]
  enum DPLEOutput {
    /// Lyapunov matrix (horzcat when const_dim, blkdiag otherwise) (Cholesky of P if pos_def) [p]
    DPLE_P,
    /// Number of arguments.
    DPLE_NUM_OUT
  };

  /// Forward declaration of internal class
  class DpleInternal;

  /**  \brief Base class for Discrete Periodic Lyapunov Equation Solvers

     @copydoc DPLE_doc
     
      \generalsection{DpleSolver}
      \pluginssection{DpleSolver}
      
       \author Joris Gillis
      \date 2014

  */
  class CASADI_CONTROL_EXPORT DpleSolver : public Function {
  public:
    /// Default constructor
    DpleSolver();

    /// Clone
    DpleSolver clone() const;

    /** \brief DpleSolver solver factory
    * \param name \pluginargument{DpleSolver}
    * \param[in] A List of sparsities of A_i
    * \param[in] V List of sparsities of V_i
    */
    DpleSolver(const std::string& name,
               const std::vector< Sparsity > & A,
               const std::vector< Sparsity > & V);

    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;

    /// Access functions of the node
    DpleInternal* operator->();

    /// Access functions of the node
    const DpleInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);
  };

} // namespace casadi

#endif // CASADI_DPLE_SOLVER_HPP
