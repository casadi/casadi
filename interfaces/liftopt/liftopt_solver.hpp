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

#ifndef LIFTOPT_SOLVER_HPP
#define LIFTOPT_SOLVER_HPP

#include <casadi/fx/nlp_solver.hpp>
#include <casadi/fx/mx_function.hpp>

namespace CasADi{
  namespace Interfaces{

    enum LOFunInputs{LO_U, LO_LAMBDA, LO_NUM_IN};
    enum LOFunOutputs{LO_OBJRES, LO_EQ, LO_INEQ, LO_OBJ, LO_LAGFCN, LO_NUM_OUT};

// Forward declaration of internal class 
class LiftoptInternal;

class LiftoptSolver : public NLPSolver{
public:

  /** \brief  Default constructor */
  LiftoptSolver();
  
  /** \brief  Create an KINSOL instance */
  explicit LiftoptSolver(const MXFunction& fcn);
  
  /** \brief  Access functions of the node */
  LiftoptInternal* operator->();

  /** \brief  Const access functions of the node */
  const LiftoptInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Node init?
  std::vector<double> &nodeInit();
  const std::vector<double> &nodeInit() const;
  
  
};

  } // namespace Interfaces
} // namespace CasADi

#endif // LIFTOPT_SOLVER_HPP
