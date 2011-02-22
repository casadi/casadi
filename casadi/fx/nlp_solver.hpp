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

#include "fx.hpp"

namespace CasADi{

/// Input arguments of an NLP Solver
enum NLPInput{
/// Decision variables initial guess
NLP_X_INIT,
/// Decision variables lower bound
NLP_LBX,
/// Decision variables upper bound
NLP_UBX,
/// Constraints lower bound
NLP_LBG,
/// Constraints upper bound
NLP_UBG,
/// Lambda multipliers initial guess
NLP_LAMBDA_INIT,
/// Lower bound multipliers initial guess
NLP_LAMBDA_LBX_INIT,
/// Upper bound multipliers initial guess
NLP_LAMBDA_UBX_INIT,
NLP_NUM_IN};

/// Outputs arguments of an NLP Solver
enum NLPOutput{
/// Decision variables for optimal solution
NLP_X_OPT,
/// Objective/cost function for optimal solution
NLP_COST,
///  Lambda multipliers function for optimal solution
NLP_LAMBDA_OPT,
///  Lower bound multipliers for optimal solution
NLP_LAMBDA_LBX,
///  Upper bound multipliers for optimal solution
NLP_LAMBDA_UBX,
NLP_NUM_OUT};

class NLPSolverInternal;

/** \brief NLPSolver

NLPSolver is an CasADi::FX mappinf from CasADi::NLPInput to CasADi::NLPOutput

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

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace CasADi

#endif // NLP_SOLVER_HPP

