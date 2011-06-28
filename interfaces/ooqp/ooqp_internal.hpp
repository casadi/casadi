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

#ifndef OOQP_INTERNAL_HPP
#define OOQP_INTERNAL_HPP

#include "casadi/fx/qp_solver_internal.hpp"

namespace CasADi{
namespace Interfaces {
  
class OOQPInternal : public QPSolverInternal {
  friend class OOQPSolver;
public:
  /** \brief  Constructor */
  explicit OOQPInternal();

  /** \brief  Clone */
  virtual OOQPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit OOQPInternal(const CRSSparsity & H, const CRSSparsity & G, const CRSSparsity & A);

  /** \brief  Destructor */
  virtual ~OOQPInternal();

  /** \brief  Initialize stage */
  virtual void init();
  
  protected:
  
};


} // namespace Interfaces
} // namespace CasADi

#endif //OOQP_INTERNAL_HPP

