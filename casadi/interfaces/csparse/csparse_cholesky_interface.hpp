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


#ifndef CASADI_CSPARSE_CHOLESKY_INTERFACE_HPP
#define CASADI_CSPARSE_CHOLESKY_INTERFACE_HPP

/** \defgroup plugin_Linsol_csparsecholesky
   * Linsol with CSparseCholesky Interface
*/

/** \pluginsection{Linsol,csparsecholesky} */

/// \cond INTERNAL
#include <cs.h>
#include "casadi/core/linsol_internal.hpp"
#include <casadi/interfaces/csparse/casadi_linsol_csparsecholesky_export.h>

namespace casadi {

  struct CASADI_LINSOL_CSPARSECHOLESKY_EXPORT CsparseCholMemory : public LinsolMemory {
    // Destructor
    ~CsparseCholMemory();

    // The transpose of linear system in form (CCS)
    cs A;

    // The symbolic factorization
    css *S;

    // The numeric factorization
    csn *L;

    // Temporary
    std::vector<double> temp;
  };

  /** \brief \pluginbrief{Linsol,csparsecholesky}
   *
   *
   @copydoc Linsol_doc
   @copydoc plugin_Linsol_csparsecholesky
   *
   */
  class CASADI_LINSOL_CSPARSECHOLESKY_EXPORT
  CSparseCholeskyInterface : public LinsolInternal {
  public:
    // Create a linear solver given a sparsity pattern and a number of right hand sides
    CSparseCholeskyInterface(const std::string& name);

    /** \brief  Create a new LinsolInternal */
    static LinsolInternal* creator(const std::string& name) {
      return new CSparseCholeskyInterface(name);
    }

    // Destructor
    ~CSparseCholeskyInterface() override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new CsparseCholMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<CsparseCholMemory*>(mem);}

    // Set sparsity pattern
    void reset(void* mem, const int* sp) const override;

    // Symbolic factorization
    void pivoting(void* mem, const double* A) const override;

    // Factorize the linear system
    void factorize(void* mem, const double* A) const override;

    // Solve the linear system
    void solve(void* mem, double* x, int nrhs, bool tr) const override;

    // Solve the system of equations <tt>Lx = b</tt>
    void solve_cholesky(void* mem, double* x, int nrhs, bool tr) const override;

    /// Obtain a symbolic Cholesky factorization
    Sparsity linsol_cholesky_sparsity(void* mem, bool tr) const override;

    /// Obtain a numeric Cholesky factorization
    DM linsol_cholesky(void* mem, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    const char* plugin_name() const override { return "csparsecholesky";}

    // Get name of the class
    std::string class_name() const override { return "CSparseCholeskyInterface";}
  };

} // namespace casadi

/// \endcond
#endif // CASADI_CSPARSE_CHOLESKY_INTERFACE_HPP
