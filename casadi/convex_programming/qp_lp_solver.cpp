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

#include "qp_lp_internal.hpp"
#include "qp_lp_solver.hpp"

using namespace std;
namespace casadi{

QPLPSolver::QPLPSolver(){
}


QPLPSolver::QPLPSolver(const LPStructure & st)  {
  assignNode(new QPLPInternal(st));
}

QPLPInternal* QPLPSolver::operator->(){
  return static_cast<QPLPInternal*>(Function::operator->());
}

const QPLPInternal* QPLPSolver::operator->() const{
  return static_cast<const QPLPInternal*>(Function::operator->());

}

bool QPLPSolver::checkNode() const{
  return dynamic_cast<const QPLPInternal*>(get());
}

QPSolver & QPLPSolver::getSolver() {
  return (*this)->qpsolver_;
}

} // namespace casadi
