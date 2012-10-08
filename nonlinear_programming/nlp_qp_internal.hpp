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

  /** \brief Internal class for NlpQPInternal
   * 
      @copydoc QPSolver_doc
   * */
class NlpQPInternal : public QPSolverInternal {
  friend class NlpQPSolver;
public:
  /** \brief  Constructor */
  explicit NlpQPInternal();

  /** \brief  Clone */
  virtual NlpQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit NlpQPInternal(const CRSSparsity& H, const CRSSparsity &A);

  /** \brief  Destructor */
  virtual ~NlpQPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate(int nfdir, int nadir);
  
  protected:
    NLPSolver nlpsolver_;
    
    // an MX that represents H, but dependant on a common MX.
    MX H_;
    // an MX that represents G, but dependant on a common MX.
    MX G_;
    // an MX that represents A, but dependant on a common MX.
    MX A_;
};

} // namespace CasADi

#endif //NLP_QP_INTERNAL_HPP

