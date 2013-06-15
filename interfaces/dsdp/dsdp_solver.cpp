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

#include "dsdp_internal.hpp"
#include "dsdp_solver.hpp"

using namespace std;
namespace CasADi{

DSDPSolver::DSDPSolver(){ 
}


DSDPSolver::DSDPSolver(const SDPStructure &st)  {
  assignNode(new DSDPInternal(st));
}

DSDPInternal* DSDPSolver::operator->(){
  return (DSDPInternal*)(FX::operator->());
}

const DSDPInternal* DSDPSolver::operator->() const{
  return (const DSDPInternal*)(FX::operator->());

}

bool DSDPSolver::checkNode() const{
  return dynamic_cast<const DSDPInternal*>(get());
}

} // namespace CasADi
