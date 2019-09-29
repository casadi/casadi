/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *    Copyright (C) 2019 Jorn Baayen, KISTERS AG
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


#ifndef CASADI_LINSOL_TRIDIAG_HPP
#define CASADI_LINSOL_TRIDIAG_HPP

/** \defgroup plugin_Linsol_tridiag
  * Linear solver for tridiagonal matrices
*/

/** \pluginsection{Linsol,tridiag} */

/// \cond INTERNAL
#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_tridiag_export.h>

namespace casadi {
  struct CASADI_LINSOL_TRIDIAG_EXPORT LinsolTridiagMemory : public LinsolMemory {
    bool have_c, have_ctr;
    std::vector<double> c, ctr, d;
  };

  /** \brief \pluginbrief{LinsolInternal,tridiag}
   * @copydoc LinsolInternal_doc
   * @copydoc plugin_LinsolInternal_tridiag
   */
  class CASADI_LINSOL_TRIDIAG_EXPORT LinsolTridiag : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    LinsolTridiag(const std::string& name, const Sparsity& sp);

    /** \brief  Create a new LinsolInternal */
    static LinsolInternal* creator(const std::string& name, const Sparsity& sp) {
      return new LinsolTridiag(name, sp);
    }

    // Destructor
    ~LinsolTridiag() override;

    // Initialize the solver
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LinsolTridiagMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LinsolTridiagMemory*>(mem);}

    // Symbolic factorization
    int nfact(void* mem, const double* A) const override;

    // Factorize the linear system
    int sfact(void* mem, const double* A) const override;

    // Solve the linear system
    int solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const override;

    /// Generate C code
    void generate(CodeGenerator& g, const std::string& A, const std::string& x,
                  casadi_int nrhs, bool tr) const override;

    // Get name of the plugin
    const char* plugin_name() const override { return "tridiag";}

    // Get name of the class
    std::string class_name() const override { return "LinsolTridiag";}

    /// A documentation string
    static const std::string meta_doc;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_LINSOL_TRIDIAG_HPP
