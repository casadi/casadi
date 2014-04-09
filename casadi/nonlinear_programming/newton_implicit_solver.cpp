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

#include "newton_implicit_internal.hpp"
#include "newton_implicit_solver.hpp"

using namespace std;
namespace casadi{

  NewtonImplicitSolver::NewtonImplicitSolver(){
  }

  NewtonImplicitSolver::NewtonImplicitSolver(const Function& f, const Function& jac, const LinearSolver& linsol)  {
    assignNode(new NewtonImplicitInternal(f,jac,linsol));
  }

  NewtonImplicitInternal* NewtonImplicitSolver::operator->(){
    return static_cast<NewtonImplicitInternal*>(Function::operator->());
  }

  const NewtonImplicitInternal* NewtonImplicitSolver::operator->() const{
    return static_cast<const NewtonImplicitInternal*>(Function::operator->());
  }

  bool NewtonImplicitSolver::checkNode() const{
    return dynamic_cast<const NewtonImplicitInternal*>(get());
  }

} // namespace casadi
