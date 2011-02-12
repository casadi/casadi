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

#ifndef KINSOL_SOLVER_HPP
#define KINSOL_SOLVER_HPP

#include "casadi/fx/implicit_function.hpp"

namespace CasADi{
namespace Sundials{
  
// Forward declaration of internal class 
class KinsolInternal;

class KinsolSolver : public ImplicitFunction{
public:

  /** \brief  Default constructor */
  KinsolSolver();
  
  /** \brief  Create an KINSOL instance */
  explicit KinsolSolver(const FX& f, int nrhs=1);
  
  /** \brief  Access functions of the node */
  KinsolInternal* operator->();

  /** \brief  Const access functions of the node */
  const KinsolInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /** \brief Generate a linear solver for the sensitivity equations */
  KinsolSolver jac(int iind=0, int oind=0);

};


} // namespace Sundials
} // namespace CasADi

#endif //KINSOL_SOLVER_HPP

