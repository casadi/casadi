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

#include "nlp_implicit_internal.hpp"
#include "nlp_implicit_solver.hpp"

using namespace std;
namespace casadi{

  NLPImplicitSolver::NLPImplicitSolver(){
  }

  NLPImplicitSolver::NLPImplicitSolver(const Function& f, const Function& jac, const LinearSolver& linsol)  {
    assignNode(new NLPImplicitInternal(f,jac,linsol));
  }

  NLPImplicitInternal* NLPImplicitSolver::operator->(){
    return static_cast<NLPImplicitInternal*>(Function::operator->());
  }

  const NLPImplicitInternal* NLPImplicitSolver::operator->() const{
    return static_cast<const NLPImplicitInternal*>(Function::operator->());
  }

  bool NLPImplicitSolver::checkNode() const{
    return dynamic_cast<const NLPImplicitInternal*>(get());
  }

  NLPSolver& NLPImplicitSolver::getNLPSolver(){
    casadi_assert(checkNode());
    return (*this)->nlp_solver_;
  }

} // namespace casadi
