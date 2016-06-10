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


#ifndef CASADI_LINSOL_HPP
#define CASADI_LINSOL_HPP

#include "function.hpp"

namespace casadi {

  // Forward declaration of internal class
  class LinsolInternal;

  /** \brief Linear solver
    * Create a solver for linear systems of equations
    * Solves the linear system A*X = B or A^T*X = B for X
    * with A square and non-singular
    *
    *  If A is structurally singular, an error will be thrown during init.
    *  If A is numerically singular, the prepare step will fail.

      \generalsection{Linsol}
      \pluginssection{Linsol}

      \author Joel Andersson
      \date 2011-2016
  */
  class CASADI_EXPORT Linsol : public SharedObject {
  public:

    /// Default constructor
    Linsol();

    /// Importer factory
    explicit Linsol(const std::string& name, const std::string& solver,
                    const Dict& opts=Dict());

    /// Access functions of the node
    LinsolInternal* operator->();
    const LinsolInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    /// Solve numerically
    DM solve(const DM& A, const DM& B, bool tr=false);

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool tr=false);

#ifndef SWIG
    // Set sparsity pattern
    void reset(const int* sp) const;

    // Select pivots
    void pivoting(const double* A) const;

    // Factorize linear system of equations
    void factorize(const double* A) const;

    // Solve factorized linear system of equations
    void solve(double* x, int nrhs=1, bool tr=false) const;

    /** \brief Solve the system of equations <tt>Lx = b</tt>
        Only when a Cholesky factorization is available
    */
    void solve_cholesky(double* x, int nrhs, bool tr) const;
#endif // SWIG

    /** \brief Obtain a symbolic Cholesky factorization
        Only for Cholesky solvers
    */
    Sparsity cholesky_sparsity(bool tr=false) const;

    /** \brief Obtain a numeric Cholesky factorization
        Only for Cholesky solvers
     */
    DM cholesky(bool tr=false) const;

    /** \brief Number of negative eigenvalues
      * Not available for all solvers
      */
    int neig() const;

    /** \brief Matrix rank
      * Not available for all solvers
      */
    int rank() const;
  };


  /// [DEPRECATED] Create a linear solver (legacy syntax, use Linsol constructor instead)
  CASADI_EXPORT Function linsol_new(const std::string& name, const std::string& solver,
                                const Sparsity& sp, int nrhs, const Dict& opts=Dict());

  /// [DEPRECATED] Legacy - alias for Linsol::has_plugin
  inline bool has_linsol(const std::string& name) { return Linsol::has_plugin(name);}

  /// [DEPRECATED] Legacy - alias for Linsol::load_plugin
  inline void load_linsol(const std::string& name) { Linsol::load_plugin(name);}

  /// [DEPRECATED] Legacy - alias for Linsol::doc
  inline std::string doc_linsol(const std::string& name) { return Linsol::doc(name);}

} // namespace casadi

#endif // CASADI_LINSOL_HPP
