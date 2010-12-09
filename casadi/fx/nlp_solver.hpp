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

#ifndef NLP_SOLVER_HPP
#define NLP_SOLVER_HPP

#include "../expression_tools.hpp"
#include "fx.hpp"

namespace CasADi{

  /// Inputs of an NLP Solver
  enum NLPInput{NLP_X_INIT,NLP_LBX,NLP_UBX,NLP_LBG,NLP_UBG,NLP_LAMBDA_INIT,NLP_NUM_IN};

  /// Outputs of an NLP Solver
  enum NLPOutput{NLP_X_OPT,NLP_COST,NLP_LAMBDA_OPT,NLP_LAMBDA_LBX,NLP_LAMBDA_UBX,NLP_NUM_OUT};

class NLPSolverInternal;

/** \brief NLPSolver
  \author Joel Andersson 
  \date 2010
*/
class NLPSolver : public FX{
  public:

  /// Default constructor
  NLPSolver();

  /// Access functions of the node
  NLPSolverInternal* operator->();
  const NLPSolverInternal* operator->() const;

  /// Assert that the node is pointing to the right type of object
  void assertNode() const;
};

} // namespace CasADi

#endif // NLP_SOLVER_HPP

