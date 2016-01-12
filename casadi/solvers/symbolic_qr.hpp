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

#include "casadi/core/function/linsol.hpp"
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

  /** \brief \pluginbrief{Linsol,symbolicqr}

      @copydoc Linsol_doc
      @copydoc plugin_Linsol_symbolicqr
      \author Joel Andersson
      \date 2013
  */
  class CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQr : public Linsol {
  public:
    // Constructor
    SymbolicQr(const std::string& name, const Sparsity& sparsity, int nrhs);

    // Destructor
    virtual ~SymbolicQr();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "symbolicqr";}

    /** \brief  Create a new Linsol */
    static Linsol* creator(const std::string& name, const Sparsity& sp, int nrhs) {
      return new SymbolicQr(name, sp, nrhs);
    }

    // Initialize
    virtual void init(const Dict& opts);

    /** \brief Allocate memory block */
    virtual Memory* memory() const;

    // Setup memory block
    virtual void set_temp(Memory& mem, const double** arg, double** res,
                          int* iw, double* w) const;

    // Factorize the linear system
    virtual void linsol_factorize(Memory& mem, const double* A) const;

    // Solve the linear system
    virtual void linsol_solve(Memory& mem, double* x, int nrhs, bool tr) const;

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(CodeGenerator& g) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(CodeGenerator& g) const;

    /** \brief Evaluate symbolically (SX) */
    virtual void linsol_eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem,
                               bool tr, int nrhs);

    // Factorization function
    Function fact_fcn_;

    // Solve function
    Function solv_fcn_N_, solv_fcn_T_;

    /// A documentation string
    static const std::string meta_doc;
  };

  /** \brief Memory for SymbolicQR  */
  struct CASADI_LINSOL_SYMBOLICQR_EXPORT SymbolicQrMemory : public WorkMemory {
    // Destructor
    virtual ~SymbolicQrMemory() {}

    // Storage for QR factorization
    std::vector<double> q, r;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SYMBOLIC_QR_HPP

