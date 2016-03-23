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


#ifndef CASADI_LINSOL_IMPL_HPP
#define CASADI_LINSOL_IMPL_HPP

#include "linsol.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Internal class
      @copydoc Linsol_doc
  */
  class CASADI_EXPORT Linsol : public FunctionInternal, public PluginInterface<Linsol> {
  public:
    /// Constructor
    Linsol(const std::string& name, const Sparsity& sparsity, int nrhs);

    /// Destructor
    virtual ~Linsol();

    ///@{
    /** \brief Number of function inputs and outputs */
    virtual size_t get_n_in() { return LINSOL_NUM_IN;}
    virtual size_t get_n_out() { return LINSOL_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs */
    virtual Sparsity get_sparsity_in(int i);
    virtual Sparsity get_sparsity_out(int i);
    /// @}

    ///@{
    /** \brief Names of function input and outputs */
    virtual std::string get_name_in(int i) { return linsol_in(i);}
    virtual std::string get_name_out(int i) { return linsol_out(i);}
    /// @}

    /// Initialize
    virtual void init(const Dict& opts);

    /// Solve the system of equations
    virtual void eval(void* mem, const double** arg, double** res, int* iw, double* w) const;

    /// Create a solve node
    using FunctionInternal::linsol_solve;
    virtual MX linsol_solve(const MX& A, const MX& B, bool tr);

    /// Evaluate SX, possibly transposed
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                               bool tr, int nrhs);

    /// Evaluate SX
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) {
      linsol_eval_sx(arg, res, iw, w, 0, false, size2_out(LINSOL_X));
    }

    /** \brief Calculate forward mode directional derivatives */
    virtual void linsol_forward(const std::vector<MX>& arg, const std::vector<MX>& res,
                                const std::vector<std::vector<MX> >& fseed,
                                std::vector<std::vector<MX> >& fsens, bool tr);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void linsol_reverse(const std::vector<MX>& arg, const std::vector<MX>& res,
                                const std::vector<std::vector<MX> >& aseed,
                                std::vector<std::vector<MX> >& asens, bool tr);

    /** \brief  Propagate sparsity forward */
    virtual void linsol_spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
                              bool tr, int nrhs);

    /** \brief  Propagate sparsity backwards */
    virtual void linsol_spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem,
                              bool tr, int nrhs);

    ///@{
    /// Propagate sparsity through a linear solve
    virtual void linsol_spsolve(bvec_t* X, const bvec_t* B, bool tr) const;
    virtual void linsol_spsolve(DM& X, const DM& B, bool tr) const;
    ///@}

    /// Dulmage-Mendelsohn decomposition
    std::vector<int> rowperm_, colperm_, rowblock_, colblock_;

    /// Get sparsity pattern
    int nrow() const { return size1_in(LINSOL_A);}
    int ncol() const { return size2_in(LINSOL_A);}
    int nnz() const { return nnz_in(LINSOL_A);}
    const int* row() const { return sparsity_in(LINSOL_A).row();}
    const int* colind() const { return sparsity_in(LINSOL_A).colind();}

    // Creator function for internal class
    typedef Linsol* (*Creator)(const std::string& name, const Sparsity& sp, int nrhs);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

    /// Infix
    static const std::string infix_;

    // Get name of the plugin
    virtual const char* plugin_name() const { return "none";}

    // Sparsity of the linear system
    Sparsity sparsity_;

    // Number of equations
    int neq_;

    // Number of right-hand-sides
    int nrhs_;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_LINSOL_IMPL_HPP
