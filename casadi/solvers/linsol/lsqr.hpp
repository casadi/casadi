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


#ifndef CASADI_LSQR_HPP
#define CASADI_LSQR_HPP

#include "casadi/core/linsol_internal.hpp"
#include <casadi/solvers/linsol/casadi_linsol_lsqr_export.h>

/** \defgroup plugin_Linsol_symbolicqr

    Linear solver for sparse least-squares problems
    Inspired from https://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/lsqr.py#L96

*/

/** \pluginsection{Linsol,symbolicqr} */

/// \cond INTERNAL

namespace casadi {

  /** \brief Memory for SymbolicQR  */
  struct CASADI_LINSOL_LSQR_EXPORT LsqrMemory : public LinsolMemory {
    // Work vectors
    std::vector<const double*> arg;
    std::vector<double> w;

    std::vector<double> A;
  };

  /** \brief \pluginbrief{Linsol,symbolicqr}

      @copydoc Linsol_doc
      @copydoc plugin_Linsol_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINSOL_LSQR_EXPORT Lsqr : public LinsolInternal {
  public:
    // Constructor
    Lsqr(const std::string& name);

    // Destructor
    ~Lsqr() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "symbolicqr";}

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name) {
      return new Lsqr(name);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    // Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_memory() const override { return new LsqrMemory();}

    /** \brief Free memory block */
    void free_memory(void *mem) const override { delete static_cast<LsqrMemory*>(mem);}

    /** \brief Initalize memory block */
    void init_memory(void* mem) const override;

    // Set sparsity pattern
    void reset(void* mem, const int* sp) const override;

    // Factorize the linear system
    void factorize(void* mem, const double* A) const override;

    // Solve the linear system
    void solve(void* mem, double* x, int nrhs, bool tr) const override;

    /// A documentation string
    static const std::string meta_doc;

    // Generated function options
    Dict fopts_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_LSQR_HPP
