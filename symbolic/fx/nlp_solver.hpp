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


/** \defgroup NLPSolver_doc

  Solves the following nonlinear optimization problem:
  \verbatim
  min          F(x,p)
   x
  
  subject to
              LBG <= G(x,p) <= UBG
              LBX <= x    <= UBX
              
      n: number of decision variables (x)
      m: number of constraints (A)
  \endverbatim

*/

namespace CasADi{
  
/// Input arguments of an NLP Solver [nlpsolverIn]
enum NLPInput{
/// Decision variables initial guess (nx x 1)  [x_init]
NLP_X_INIT,
/// Decision variables lower bound (nx x 1), default -inf [lbx]
NLP_LBX,
/// Decision variables upper bound (nx x 1), default +inf [ubx]
NLP_UBX,
/// Constraints lower bound (ng x 1), default -inf [lbg]
NLP_LBG,
/// Constraints upper bound (ng x 1), default +inf [ubg]
NLP_UBG,
/// Lagrange multipliers associated with G, initial guess (ng x 1) [lambda_init]
NLP_LAMBDA_INIT,
/// Parameters on which the objective and constraints might depend (np x 1) [p]
NLP_P,
NLP_NUM_IN};

/// Output arguments of an NLP Solver [nlpsolverOut]
enum NLPOutput{
/// Decision variables for optimal solution (nx x 1) [x_opt]
NLP_X_OPT,
/// Objective/cost function for optimal solution (1 x 1) [cost]
NLP_COST,
/// Lagrange multipliers associated with G at the solution (ng x 1) [lambda_g]
NLP_LAMBDA_G,
/// Lagrange multipliers associated with bounds on X at the solution (nx x 1) [lambda_x]
NLP_LAMBDA_X, 
/// Lagrange multipliers associated with the parameters (np x 1) [lambda_p]
NLP_LAMBDA_P, 
/// The constraints evaluated at the optimal solution (ng x 1) [g]
NLP_G,
NLP_NUM_OUT};

class NLPSolverInternal;

/** \brief NLPSolver

  @copydoc NLPSolver_doc

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
  
  /// Prints out a human readable report about possible constraint violations, after solving 
  void reportConstraints(std::ostream &stream=std::cout);

  std::string getReportConstraints() { std::stringstream s; reportConstraints(s); return s.str(); }
  
  /// Set options that make the NLP solver more suitable for solving QPs
  void setQPOptions();
  
  /// Access the objective function F
  FX getF() const;
  
  /// Access the objective function G
  FX getG() const;

  /// Access the hessian of the Lagrangian function H
  FX getH() const;
  
  /// Access the jacobian of the constraint function J
  FX getJ() const;
    
};

} // namespace CasADi

#endif // NLP_SOLVER_HPP

