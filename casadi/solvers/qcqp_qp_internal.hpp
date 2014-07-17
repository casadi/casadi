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

#ifndef CASADI_QCQP_QP_INTERNAL_HPP
#define CASADI_QCQP_QP_INTERNAL_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include "casadi/core/function/qcqp_solver.hpp"
#include <casadi/solvers/casadi_qpsolver_qcqp_export.h>

/// \cond INTERNAL
namespace casadi {

  /** \brief Use a QCQP solver to solve q QP

      @copydoc QpSolver_doc

      \author Joris Gillis
      \date 2013
  */
  class CASADI_QPSOLVER_QCQP_EXPORT QCQPQPInternal : public QpSolverInternal {
  public:

    /** \brief  Create a new Solver */
    explicit QCQPQPInternal(const std::vector<Sparsity> &st);

    /** \brief  Destructor */
    virtual ~QCQPQPInternal();

    /** \brief  Clone */
    virtual QCQPQPInternal* clone() const;

    /** \brief  Create a new QP Solver */
    static QpSolverInternal* creator(const QPStructure& st)
    { return new QCQPQPInternal(st);}

    /** \brief  Initialize */
    virtual void init();

    virtual void evaluate();

  protected:
    QcqpSolver qcqpsolver_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_QCQP_QP_INTERNAL_HPP
