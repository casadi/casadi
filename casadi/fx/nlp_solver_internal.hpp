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

#ifndef NLP_SOLVER_INTERNAL_HPP
#define NLP_SOLVER_INTERNAL_HPP

#include "nlp_solver.hpp"
#include "fx_internal.hpp"

namespace CasADi{
    
/** \brief NLP solver storage class

  @copydoc NLPSolver_doc
  \author Joel Andersson 
  \date 2010
*/
class NLPSolverInternal : public FXInternal{

public:
  explicit NLPSolverInternal();
  virtual ~NLPSolverInternal() = 0;

  virtual void init();

protected:
  int n_,m_;
  
  int pn_,pm_; // (n x m) parameterspace
  
};

} // namespace CasADi

#endif //NLP_SOLVER_INTERNAL_HPP
