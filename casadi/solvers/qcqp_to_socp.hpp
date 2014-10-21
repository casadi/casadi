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


#ifndef CASADI_QCQP_TO_SOCP_HPP
#define CASADI_QCQP_TO_SOCP_HPP

#include "casadi/core/function/qcqp_solver_internal.hpp"
#include "casadi/core/function/socp_solver_internal.hpp"
#include "casadi/core/function/linear_solver.hpp"
#include <casadi/solvers/casadi_qcqpsolver_socp_export.h>

/** \defgroup plugin_QcqpSolver_socp

   * Solve a QCQP with an SocpSolver
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

*/

/** \pluginsection{QcqpSolver,socp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{QcqpSolver,socp}

   @copydoc QcqpSolver_doc
   @copydoc plugin_QcqpSolver_socp

   \author Joris Gillis
   \date 2013
  */
  class CASADI_QCQPSOLVER_SOCP_EXPORT QcqpToSocp : public QcqpSolverInternal,
    public Adaptor<QcqpToSocp, SocpSolverInternal> {
  public:

    /** \brief  Create a new Solver */
    explicit QcqpToSocp(const std::vector<Sparsity> &st);

    /** \brief  Destructor */
    virtual ~QcqpToSocp();

    /** \brief  Clone */
    virtual QcqpToSocp* clone() const;

    /** \brief  Create a new QP Solver */
    static QcqpSolverInternal* creator(const QCQPStructure& st)
    { return new QcqpToSocp(st);}

    /** \brief  Initialize */
    virtual void init();

    virtual void evaluate();

    /// A documentation string
    static const std::string meta_doc;

    /// Solve with
    SocpSolver solver_;

  protected:
    std::vector<LinearSolver> cholesky_;
  };
  /// \endcond
} // namespace casadi

#endif // CASADI_QCQP_TO_SOCP_HPP
