/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "printable.hpp"

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

      \identifier{1kh} */
  class CASADI_EXPORT Linsol
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Linsol>) {
  public:
    /** \brief Get type name

        \identifier{1ki} */
    static std::string type_name() {return "Linsol";}

    /// Default constructor
    Linsol();

    /// Constructor
    explicit Linsol(const std::string& name, const std::string& solver,
                    const Sparsity& sp, const Dict& opts=Dict());

    /// Access functions of the node
    LinsolInternal* operator->();
    const LinsolInternal* operator->() const;

#ifndef SWIG
    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);
#endif //SWIG

    /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{1kj} */
    static Linsol create(LinsolInternal* node);
#endif // SWIG
    /// \endcond

    /// Check if a plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    /// Get linear system sparsity
    const Sparsity& sparsity() const;

    /// Symbolic factorization of the linear system, e.g. selecting pivots
    void sfact(const DM& A) const;

    /// Numeric factorization of the linear system
    void nfact(const DM& A) const;

    ///@{
    /// Solve linear system of equations
    DM solve(const DM& A, const DM& B, bool tr=false) const;
    MX solve(const MX& A, const MX& B, bool tr=false) const;
    ///@}

    /** \brief Number of negative eigenvalues

      * Not available for all solvers

        \identifier{1kk} */
    casadi_int neig(const DM& A) const;

    /** \brief Matrix rank

      * Not available for all solvers

        \identifier{1kl} */
    casadi_int rank(const DM& A) const;

   /// Get all statistics obtained at the end of the last evaluate call
    Dict stats(int mem=1) const;

    #ifndef SWIG
    ///@{
    /// Low-level API
    int sfact(const double* A, int mem=0) const;
    int nfact(const double* A, int mem=0) const;
    int solve(const double* A, double* x, casadi_int nrhs=1, bool tr=false, int mem=0) const;
    casadi_int neig(const double* A, int mem=0) const;
    casadi_int rank(const double* A, int mem=0) const;
    ///@}

    /// Checkout a memory object
    casadi_int checkout() const;

    /// Release a memory object
    void release(int mem) const;

    /** \brief Serialize an object

        \identifier{1km} */
    void serialize(SerializingStream &s) const;

    /** \brief Deserialize with type disambiguation

        \identifier{1kn} */
    static Linsol deserialize(DeserializingStream& s);

    #endif // SWIG
  };

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_linsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_linsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_linsol(const std::string& name);

} // namespace casadi

#endif // CASADI_LINSOL_HPP
