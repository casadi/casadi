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


#ifndef CASADI_QP_TO_NLP_HPP
#define CASADI_QP_TO_NLP_HPP

#include "casadi/core/function/qp_solver_internal.hpp"
#include "casadi/core/function/nlp_solver_internal.hpp"

#include <casadi/solvers/casadi_qpsolver_nlp_export.h>


/** \defgroup plugin_QpSolver_nlp
   Solve QPs using an NlpSolver
*/

/** \pluginsection{QpSolver,nlp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{QpSolver,nlp}

   @copydoc QpSolver_doc
   @copydoc plugin_QpSolver_nlp

   \author Joris Gillis
   \date 2011
  */
class CASADI_QPSOLVER_NLP_EXPORT QpToNlp : public QpSolverInternal,
  public Adaptor<QpToNlp, NlpSolverInternal> {
public:
  /** \brief  Constructor */
  explicit QpToNlp();

  /** \brief  Clone */
  virtual QpToNlp* clone() const;

  /** \brief  Create a new QP Solver */
  static QpSolverInternal* creator(const QPStructure& st)
  { return new QpToNlp(st);}

  /** \brief  Create a new Solver */
  explicit QpToNlp(const std::vector<Sparsity> &st);

  /** \brief  Destructor */
  virtual ~QpToNlp();

  /** \brief  Initialize */
  virtual void init();

  virtual void evaluate();

  /// A documentation string
  static const std::string meta_doc;

  /// Solve with
  NlpSolver solver_;
};

} // namespace casadi
/// \endcond
#endif // CASADI_QP_TO_NLP_HPP
