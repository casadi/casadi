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

#include "socp_qcqp_internal.hpp"
#include "socp_qcqp_solver.hpp"

using namespace std;
namespace CasADi{

SOCPQCQPSolver::SOCPQCQPSolver(){ 
}


SOCPQCQPSolver::SOCPQCQPSolver(const QCQPStructure & st)  {
  assignNode(new SOCPQCQPInternal(st));
}

SOCPQCQPInternal* SOCPQCQPSolver::operator->(){
  return (SOCPQCQPInternal*)(Function::operator->());
}

const SOCPQCQPInternal* SOCPQCQPSolver::operator->() const{
  return (const SOCPQCQPInternal*)(Function::operator->());

}

bool SOCPQCQPSolver::checkNode() const{
  return dynamic_cast<const SOCPQCQPInternal*>(get());
}

SOCPSolver & SOCPQCQPSolver::getSolver() {
  return (*this)->socpsolver_;
}

} // namespace CasADi
