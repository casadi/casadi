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

#ifndef CASADI_SOCP_QCQP_INTERNAL_HPP
#define CASADI_SOCP_QCQP_INTERNAL_HPP

#include "casadi/core/function/qcqp_solver_internal.hpp"
#include "casadi/core/function/socp_solver.hpp"
#include "casadi/core/function/linear_solver.hpp"
#include <casadi/solvers/casadi_qcqpsolver_socp_export.h>

/// \cond INTERNAL
namespace casadi {

  /** \brief SOCP QCQP Solver for quadratic programming
   *
   *  Note: this implementation relies on Cholesky decomposition:
   *        <tt>Chol(H) = L ->  H = LL'</tt> with L lower triangular
   *   This requires Pi, H to be positive definite. Positive semi-definite is not sufficient.
   *    Notably, H==0  will not work.
   *
   *  A better implementation would rely on matrix square root, but we need
   *  singular value decomposition to implement that.
   *
   *
   * This implementation makes use of the epigraph reformulation:
   * \verbatim 
   *  min f(x)
   *    x
   *
   *   min  t
   *    x, t  f(x) <= t
   * \endverbatim
   *
   *  This implementation makes use of the following identity:
   * \verbatim
   *  || Gx+h||_2 <= e'x + f
   *
   *  x'(G'G - ee')x + (2 h'G - 2 f e') x + h'h - f <= 0
   * \endverbatim
   *    where we put e = [0 0 ... 1]  for the quadratic constraint
   *    arising from the epigraph reformulation and e==0 for all other 
   *    quadratic constraints.

   @copydoc QcqpSolver_doc

   \author Joris Gillis
   \date 2013
  */
  class CASADI_QCQPSOLVER_SOCP_EXPORT SOCPQCQPInternal : public QcqpSolverInternal {
  public:

    /** \brief  Create a new Solver */
    explicit SOCPQCQPInternal(const std::vector<Sparsity> &st);

    /** \brief  Destructor */
    virtual ~SOCPQCQPInternal();

    /** \brief  Clone */
    virtual SOCPQCQPInternal* clone() const;

    /** \brief  Create a new QP Solver */
    static QcqpSolverInternal* creator(const QCQPStructure& st)
    { return new SOCPQCQPInternal(st);}

    /** \brief  Initialize */
    virtual void init();

    virtual void evaluate();

  protected:
    SocpSolver socpsolver_;
    std::vector<LinearSolver> cholesky_;
  };
  /// \endcond
} // namespace casadi

#endif // CASADI_SOCP_QCQP_INTERNAL_HPP
