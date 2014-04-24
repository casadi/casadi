/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SYMBOLIC_QR_INTERNAL_HPP
#define SYMBOLIC_QR_INTERNAL_HPP

#include "symbolic_qr.hpp"
#include "linear_solver_internal.hpp"

/// \cond INTERNAL

namespace casadi {

  class CASADI_CORE_EXPORT SymbolicQRInternal : public LinearSolverInternal {
  public:
    // Constructor
    SymbolicQRInternal(const Sparsity& sparsity, int nrhs);

    // Destructor
    virtual ~SymbolicQRInternal();

    /** \brief  Clone */
    virtual SymbolicQRInternal* clone() const { return new SymbolicQRInternal(*this);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    // Initialize
    virtual void init();

    // Prepare the factorization
    virtual void prepare();

    // Solve the system of equations
    virtual void solve(double* x, int nrhs, bool transpose);

    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(std::ostream &stream, const std::string& type,
                                      CodeGenerator& gen) const;

    /** \brief Generate code for the body of the C function */
    virtual void generateBody(std::ostream &stream, const std::string& type,
                              CodeGenerator& gen) const;

    /** \brief Evaluate symbolically (SX) */
    virtual void evaluateSXGen(const SXPtrV& input, SXPtrV& output, bool tr);

    // Factorization function
    Function fact_fcn_;

    // Solve function
    Function solv_fcn_N_, solv_fcn_T_;

    // Storage for QR factorization
    DMatrix Q_, R_;
  };

} // namespace casadi

/// \endcond
#endif //SYMBOLIC_QR_INTERNAL_HPP

