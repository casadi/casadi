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


#ifndef CASADI_LR_DPLE_SOLVER_HPP
#define CASADI_LR_DPLE_SOLVER_HPP

#include "function.hpp"

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
#ifndef SWIG

  /// Input arguments of a \e dple solver [lrdpleIn]
  enum LR_DPLEInput {
    /// A matrices (horzcat when const_dim, diagcat otherwise) [a]
    LR_DPLE_A,
    /// V matrices (horzcat when const_dim, diagcat otherwise) [v]
    LR_DPLE_V,
    /// C matrix [c]
    LR_DPLE_C,
    /// H matrix: horizontal stack of all Hi [h]
    LR_DPLE_H,
    LR_DPLE_NUM_IN
  };

  /// Output arguments of a \e dple solver [lrdpleOut]
  enum LR_DPLEOutput {
    /// Lyapunov matrix (horzcat when const_dim, diagcat otherwise) (Cholesky of P if pos_def) [y]
    LR_DPLE_Y,
    /// Number of arguments.
    LR_DPLE_NUM_OUT
  };
#endif // SWIG

  /// Forward declaration of internal class
  class LrDpleInternal;

  /**  \brief Base class for Discrete Periodic Lyapunov Equation Solvers

     @copydoc DPLE_doc

      \generalsection{LrDpleSolver}
      \pluginssection{LrDpleSolver}

       \author Joris Gillis
      \date 2014

  */
  class CASADI_EXPORT LrDpleSolver : public Function {
  public:
    /// Default constructor
    LrDpleSolver();

    /// Clone
    LrDpleSolver clone() const;

    /** \brief Constructor (new syntax, includes initialization)
     * \param solver \pluginargument{LrDpleSolver}
     * \param st \structargument{LrDple}
     */
    LrDpleSolver(const std::string& name, const std::string& solver,
                 const std::map<std::string, std::vector<Sparsity> >& st,
                 const Dict& opts=Dict());

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Constructor (no initialization)
     * \param solver \pluginargument{LrDpleSolver}
     * \param st \structargument{LrDple}
     */
    LrDpleSolver(const std::string& solver,
                 const std::map<std::string, std::vector<Sparsity> >& st);
#endif // WITH_DEPRECATED_FEATURES

    /// Print solver statistics
    void printStats(std::ostream &stream=casadi::userOut()) const;

    /// Access functions of the node
    LrDpleInternal* operator->();

    /// Access functions of the node
    const LrDpleInternal* operator->() const;

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

#endif // CASADI_LR_DPLE_SOLVER_HPP
