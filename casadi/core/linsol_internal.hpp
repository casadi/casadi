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
    // Current state of factorization
    bool is_pivoted, is_factorized;

    // Constructor
    LinsolMemory() : is_pivoted(false), is_factorized(false) {}
  };

  /** Internal class
      @copydoc Linsol_doc
  */
  class CASADI_EXPORT LinsolInternal
    : public ProtoFunction, public PluginInterface<LinsolInternal> {
  public:
    /// Constructor
    LinsolInternal(const std::string& name, const Sparsity& sp);

    /// Destructor
    ~LinsolInternal() override;

    /** \brief Display object */
    void disp(std::ostream& stream, bool more) const override;

    /** \brief  Print more */
    virtual void disp_more(std::ostream& stream) const {}

    /// Initialize
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_mem() const override { return new LinsolMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override { delete static_cast<LinsolMemory*>(mem);}

    /// Evaluate SX, possibly transposed
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, void* mem,
                               bool tr, int nrhs) const;

    /// Solve Cholesky
    virtual void solve_cholesky(void* mem, double* x, int nrhs, bool tr) const;

    // Symbolic factorization - partial pivoting (optional)
    virtual void pivoting(void* mem, const double* A) const {}

    /// Factorize the linear system
    virtual void factorize(void* mem, const double* A) const;

    // Solve numerically
    virtual void solve(void* mem, double* x, int nrhs, bool tr) const;

    /// Sparsity pattern of the cholesky factors
    virtual Sparsity linsol_cholesky_sparsity(void* mem, bool tr) const;

    /// Get Cholesky factor
    virtual DM linsol_cholesky(void* mem, bool tr) const;

    /// Number of negative eigenvalues
    virtual int neig(void* mem) const;

    /// Matrix rank
    virtual int rank(void* mem) const;

    // Creator function for internal class
    typedef LinsolInternal* (*Creator)(const std::string& name, const Sparsity& sp);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    // Get name of the plugin
    const char* plugin_name() const override = 0;

    /// Get sparsity pattern
    int nrow() const { return sp_.size1();}
    int ncol() const { return sp_.size2();}
    const int* colind() const { return sp_.colind();}
    const int* row() const { return sp_.row();}
    int nnz() const { return sp_.nnz();}

    // Sparsity pattern of the linear system
    Sparsity sp_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_LINSOL_INTERNAL_HPP
