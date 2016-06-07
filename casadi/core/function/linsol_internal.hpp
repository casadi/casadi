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


#ifndef CASADI_LINSOL_INTERNAL_HPP
#define CASADI_LINSOL_INTERNAL_HPP

#include "linsol.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  struct CASADI_EXPORT LinsolMemory {
    // Sparsity pattern (allowed to change)
    std::vector<int> sparsity;
  };

  /** Internal class
      @copydoc Linsol_doc
  */
  class CASADI_EXPORT LinsolInternal
    : public FunctionInternal, public PluginInterface<LinsolInternal> {
  public:
    /// Constructor
    LinsolInternal(const std::string& name, const Sparsity& sparsity);

    /// Destructor
    virtual ~LinsolInternal();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return 0;}
    virtual size_t get_n_out() { return 0;}
    ///@}

    /// Initialize
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new LinsolMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<LinsolMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Solve numerically
    virtual void linsol_solve(void* mem, double* x, int nrhs, bool tr) const;

    /// Evaluate SX, possibly transposed
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                               bool tr, int nrhs);

    /// Solve Cholesky
    virtual void linsol_solveL(void* mem, double* x, int nrhs, bool tr) const;

    /// Factorize the linear system
    virtual void linsol_factorize(void* mem, const double* A) const;

    /// Sparsity pattern of the cholesky factors
    virtual Sparsity linsol_cholesky_sparsity(void* mem, bool tr) const;

    /// Get Cholesky factor
    virtual DM linsol_cholesky(void* mem, bool tr) const;

    /// Get sparsity pattern
    int nrow() const { return sparsity_.size1();}
    int ncol() const { return sparsity_.size2();}
    int nnz() const { return sparsity_.nnz();}
    const int* sparsity() const { return sparsity_;}
    const int* row() const { return sparsity_.row();}
    const int* colind() const { return sparsity_.colind();}

    // Creator function for internal class
    typedef LinsolInternal* (*Creator)(const std::string& name, const Sparsity& sp);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    // Get name of the plugin
    virtual const char* plugin_name() const = 0;

  private:
    // Sparsity of the linear system
    Sparsity sparsity_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_LINSOL_INTERNAL_HPP
