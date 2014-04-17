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

#ifndef SYMBOLIC_QR_HPP
#define SYMBOLIC_QR_HPP

#include "linear_solver.hpp"

namespace casadi{

  // Forward declaration of internal class
  class SymbolicQRInternal;

  /** \brief  LinearSolver based on QR factorization with sparsity pattern based reordering
             _without_ partial pivoting
      @copydoc LinearSolver_doc
      \author Joel Andersson
      \date 2013
  */
  class CASADI_SYMBOLIC_EXPORT SymbolicQR : public LinearSolver{
  public:

    /// Default (empty) constructor
    SymbolicQR();

    /// Create a linear solver given a sparsity pattern
    explicit SymbolicQR(const Sparsity& sp, int nrhs=1);

    /// Access functions of the node
    SymbolicQRInternal* operator->();

    /// Const access functions of the node
    const SymbolicQRInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static LinearSolver creator(const Sparsity& sp, int nrhs){ return SymbolicQR(sp,nrhs);}
#ifdef SWIG
    %nocallback;
#endif

  };

} // namespace casadi

#endif //SYMBOLIC_QR_HPP
