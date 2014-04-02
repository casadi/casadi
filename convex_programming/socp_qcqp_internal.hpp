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

#ifndef SOCP_QCQP_INTERNAL_HPP
#define SOCP_QCQP_INTERNAL_HPP

#include "symbolic/function/qcqp_solver_internal.hpp"
#include "symbolic/function/socp_solver.hpp"
#include "interfaces/csparse/csparse_cholesky.hpp"

/// \cond INTERNAL
namespace CasADi{

  /** \brief Internal class for SOCPQCQPInternal
   * 
      @copydoc QCQPSolver_doc
   * */
class SOCPQCQPInternal : public QCQPSolverInternal {
  friend class SOCPQCQPSolver;
public:

  /** \brief  Clone */
  virtual SOCPQCQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit SOCPQCQPInternal(const std::vector<Sparsity> &st);

  /** \brief  Destructor */
  virtual ~SOCPQCQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate();
  
  protected:
    SOCPSolver socpsolver_;
    std::vector<CSparseCholesky> cholesky_;
};
/// \endcond
} // namespace CasADi

#endif //SOCP_QCQP_INTERNAL_HPP

