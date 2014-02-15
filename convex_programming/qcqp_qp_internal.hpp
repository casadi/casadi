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

#ifndef QCQP_QP_INTERNAL_HPP
#define QCQP_QP_INTERNAL_HPP

#include "symbolic/fx/qp_solver_internal.hpp"
#include "symbolic/fx/qcqp_solver.hpp"

namespace CasADi{

  /** \brief Internal class for QCQPQPInternal
   * 
      @copydoc QPSolver_doc
   * */
class QCQPQPInternal : public QPSolverInternal {
  friend class QCQPQPSolver;
public:

  /** \brief  Clone */
  virtual QCQPQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit QCQPQPInternal(const std::vector<CRSSparsity> &st);

  /** \brief  Destructor */
  virtual ~QCQPQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate();
  
  protected:
    QCQPSolver qcqpsolver_;

};

} // namespace CasADi

#endif //QCQP_QP_INTERNAL_HPP

