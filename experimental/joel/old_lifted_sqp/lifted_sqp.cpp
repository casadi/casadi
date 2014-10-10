/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#include "lifted_sqp_internal.hpp"

using namespace std;

namespace casadi{

LiftedSQP::LiftedSQP(){
}
  
LiftedSQP::LiftedSQP(const Function& F, const Function& G){
  assignNode(new LiftedSQPInternal(F,G));
}

LiftedSQPInternal* LiftedSQP::operator->(){
  return (LiftedSQPInternal*)(NlpSolver::operator->());
}

const LiftedSQPInternal* LiftedSQP::operator->() const{
  return (const LiftedSQPInternal*)(NlpSolver::operator->());
}
    
bool LiftedSQP::checkNode() const{
  return dynamic_cast<const LiftedSQPInternal*>(get());
}

} // namespace casadi
