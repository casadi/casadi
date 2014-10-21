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


#ifndef CASADI_LP_TO_QP_HPP
#define CASADI_LP_TO_QP_HPP

#include "casadi/core/function/lp_internal.hpp"
#include "casadi/core/function/qp_solver_internal.hpp"

#include <casadi/solvers/casadi_lpsolver_qp_export.h>

/** \defgroup plugin_LpSolver_qp

    Solve LPs using a QpSolver
*/

/** \pluginsection{LpSolver,qp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{LpSolver,qp}

   @copydoc LpSolver_doc
   @copydoc plugin_LpSolver_qp

   \author Joris Gillis
   \date 2013
  */
class CASADI_LPSOLVER_QP_EXPORT LpToQp : public LpSolverInternal,
  public Adaptor<LpToQp, QpSolverInternal> {
public:

  /** \brief  Create a new Solver */
  explicit LpToQp(const std::vector<Sparsity> &st);

  /** \brief  Destructor */
  virtual ~LpToQp();

  /** \brief  Clone */
  virtual LpToQp* clone() const;

  /** \brief  Create a new QP Solver */
  static LpSolverInternal* creator(const LPStructure& st)
  { return new LpToQp(st);}

  /** \brief  Initialize */
  virtual void init();

  virtual void evaluate();

  /// A documentation string
  static const std::string meta_doc;

  /// Solve with
  QpSolver solver_;
};

} // namespace casadi
/// \endcond

#endif // CASADI_LP_TO_QP_HPP

