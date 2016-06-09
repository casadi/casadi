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


#ifndef CASADI_SYMBOLIC_QR_HPP
#define CASADI_SYMBOLIC_QR_HPP

#include "casadi/core/function/linsol_internal.hpp"
#include <casadi/solvers/casadi_linsol_symbolicqr_export.h>

/** \defgroup plugin_Linsol_symbolicqr

       Linsol based on QR factorization with sparsity pattern based reordering
      _without_ partial pivoting
*/

/** \pluginsection{Linsol,symbolicqr} */

/// \cond INTERNAL

namespace casadi {
  typedef SX* SXPtr;
  typedef std::vector<SXPtr> SXPtrV;

  /** \brief Memory for SymbolicQR  */
  struct CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQrMemory : public LinsolMemory {
    // Functions for factorization and (optionally transposed) solve
    Function factorize, solve, solveT;

    // Work vectors
    std::vector<const double*> arg;
    std::vector<double*> res;
    std::vector<int> iw;
    std::vector<double> w;

    // Allocate memory for a function
    void alloc(const Function& f);

    // Storage for QR factorization
    std::vector<double> q, r;
  };

  /** \brief \pluginbrief{Linsol,symbolicqr}

      @copydoc Linsol_doc
      @copydoc plugin_Linsol_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQr : public LinsolInternal {
  public:
    // Constructor
    SymbolicQr(const std::string& name);

    // Destructor
    virtual ~SymbolicQr();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "symbolicqr";}

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name) {
      return new SymbolicQr(name);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new SymbolicQrMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<SymbolicQrMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Set sparsity pattern
    virtual void reset(void* mem, const int* sp) const;

    // Factorize the linear system
    virtual void factorize(void* mem, const double* A) const;

    // Solve the linear system
    virtual void solve(void* mem, double* x, int nrhs, bool tr) const;

    /** \brief Evaluate symbolically (SX) */
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                               bool tr, int nrhs);

    /// A documentation string
    static const std::string meta_doc;

    // Generated function options
    Dict fopts_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SYMBOLIC_QR_HPP
