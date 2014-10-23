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


#ifndef CASADI_QP_TO_QCQP_HPP
#define CASADI_QP_TO_QCQP_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include "casadi/core/function/qcqp_solver_internal.hpp"
#include <casadi/solvers/casadi_qpsolver_qcqp_export.h>

/** \defgroup plugin_QpSolver_qcqp
   Solve QP using a QcqpSolver
*/

/** \pluginsection{QpSolver,qcqp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{QpSolver,qcqp}

      @copydoc QpSolver_doc
      @copydoc plugin_QpSolver_qcqp

      \author Joris Gillis
      \date 2013
  */
  class CASADI_QPSOLVER_QCQP_EXPORT QpToQcqp : public QpSolverInternal,
    public Adaptor<QpToQcqp, QcqpSolverInternal> {
  public:

    /** \brief  Create a new Solver */
    explicit QpToQcqp(const std::vector<Sparsity> &st);

    /** \brief  Destructor */
    virtual ~QpToQcqp();

    /** \brief  Clone */
    virtual QpToQcqp* clone() const;

    /** \brief  Create a new QP Solver */
    static QpSolverInternal* creator(const QPStructure& st)
    { return new QpToQcqp(st);}

    /** \brief  Initialize */
    virtual void init();

    virtual void evaluate();

    /// A documentation string
    static const std::string meta_doc;

    /// Solve with
    QcqpSolver solver_;
  };

} // namespace casadi
/// \endcond
#endif // CASADI_QP_TO_QCQP_HPP
