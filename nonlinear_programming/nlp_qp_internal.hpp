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

#ifndef NLP_QP_INTERNAL_HPP
#define NLP_QP_INTERNAL_HPP

#include "symbolic/fx/qp_solver_internal.hpp"
#include "symbolic/fx/nlp_solver.hpp"

namespace CasADi{

  /** \brief Internal class for NLPQPInternal
   * 
      @copydoc QPSolver_doc
   * */
class NLPQPInternal : public QPSolverInternal {
  friend class NLPQPSolver;
public:
  /** \brief  Constructor */
  explicit NLPQPInternal();

  /** \brief  Clone */
  virtual NLPQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit NLPQPInternal(const std::vector<Sparsity> &st);

  /** \brief  Destructor */
  virtual ~NLPQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate();
  
  protected:
    NLPSolver nlpsolver_;

};

} // namespace CasADi

#endif //NLP_QP_INTERNAL_HPP

