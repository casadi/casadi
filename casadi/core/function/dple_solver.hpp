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


#ifndef CASADI_DPLE_SOLVER_HPP
#define CASADI_DPLE_SOLVER_HPP

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

  /// Input arguments of a \e dple solver [dpleIn]
  enum DPLEInput {
    /// A matrices (horzcat when const_dim, diagcat otherwise) [a]
    DPLE_A,
    /// V matrices (horzcat when const_dim, diagcat otherwise) [v]
    DPLE_V,
    DPLE_NUM_IN
  };

  /// Output arguments of a \e dple solver [dpleOut]
  enum DPLEOutput {
    /// Lyapunov matrix (horzcat when const_dim, diagcat otherwise) (Cholesky of P if pos_def) [p]
    DPLE_P,
    /// Number of arguments.
    DPLE_NUM_OUT
  };
#endif // SWIG

  /// Forward declaration of internal class
  class DpleInternal;

  /**  \brief Base class for Discrete Periodic Lyapunov Equation Solvers

     @copydoc DPLE_doc

      \generalsection{DpleSolver}
      \pluginssection{DpleSolver}

       \author Joris Gillis
      \date 2014

  */
  class CASADI_EXPORT DpleSolver : public Function {
  public:
    /// Default constructor
    DpleSolver();

    /// Clone
    DpleSolver clone() const;

    /** \brief Constructor (new syntax, includes initialization)
     * \param solver \pluginargument{DpleSolver}
     * \param st \structargument{Dple}
     */
    DpleSolver(const std::string& name, const std::string& solver,
               const std::map<std::string, std::vector<Sparsity> >& st,
               const Dict& opts=Dict());

#ifdef WITH_DEPRECATED_FEATURES
    /** \brief [DEPRECATED] Constructor (no initialization)
     * \param solver \pluginargument{DpleSolver}
     * \param st \structargument{Dple}
     */
    DpleSolver(const std::string& solver,
               const std::map<std::string, std::vector<Sparsity> >& st);
#endif // WITH_DEPRECATED_FEATURES

    /// Print solver statistics
    void printStats(std::ostream &stream=casadi::userOut()) const;

    /// Access functions of the node
    DpleInternal* operator->();

    /// Access functions of the node
    const DpleInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /** \brief Obtain Periodic Schur Form of a set of matrices
     *
     *  Finds Z_i such that
     \verbatim
     Z_1' * A_1 * Z_2 = T_1,
     Z_2' * A_2 * Z_3 = T_2,
     ...
     Z_K' * A_K * Z_1 = T_K,
     \endverbatim
     *
     *  with T_1 in Hessenberg form (upper triangular + one band below 
     *  the diagonal) and T_2..T_K  upper diagonal
     *
     *  with <tt>Z_k Z_k' = eye(n) = Z_k' Z_k</tt>
     *
     */
    static void periodic_schur(const std::string& name,
                               const std::vector< Matrix<double> > & A,
                               std::vector< Matrix<double> > & SWIG_OUTPUT(T),
                               std::vector< Matrix<double> > & SWIG_OUTPUT(Z),
                               std::vector<double> &SWIG_OUTPUT(eig_real),
                               std::vector<double> &SWIG_OUTPUT(eig_imag),
                               double num_zero=0);
  };

} // namespace casadi

#endif // CASADI_DPLE_SOLVER_HPP
