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

#ifndef IPOPT_QP_SOLVER_HPP
#define IPOPT_QP_SOLVER_HPP

#include "casadi/fx/qp_solver.hpp"

namespace CasADi {
namespace Interfaces {
  
  
// Forward declaration of internal class 
class IpoptQPInternal;

  /** OOQP Solver for quadratic programing:

  min          x'Hx + G'x 
  
  subject to
              LBA <= Ax <= UBA
              LBX <= x  <= UBX
              
      nx: number of decision variables (x)
      nc: number of constraints (A)
      
      The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.
      
      NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
  */
class IpoptQPSolver : public QPSolver {
public:

  /** \brief  Default constructor */
  IpoptQPSolver();
  
  explicit IpoptQPSolver(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A);
  
  /** \brief  Access functions of the node */
  IpoptQPInternal* operator->();
  const IpoptQPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
};


} // namespace Interfaces
} // namespace CasADi

#endif //IPOPT_QP_SOLVER_HPP

