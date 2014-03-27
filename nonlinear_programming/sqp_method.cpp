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

#include "sqp_internal.hpp"

using namespace std;

namespace CasADi{

  SQPMethod::SQPMethod(){
  }
  
  SQPMethod::SQPMethod(const Function& F, const Function& G){
    assignNode(new SQPInternal(joinFG(F,G)));
  }

  SQPMethod::SQPMethod(const Function& nlp){
    assignNode(new SQPInternal(nlp));
  }

  SQPInternal* SQPMethod::operator->(){
    return static_cast<SQPInternal*>(NLPSolver::operator->());
  }

  const SQPInternal* SQPMethod::operator->() const{
    return static_cast<const SQPInternal*>(NLPSolver::operator->());
  }
    
  bool SQPMethod::checkNode() const{
    return dynamic_cast<const SQPInternal*>(get())!=0;
  }

  const QPSolver SQPMethod::getQPSolver() const {
    return (*this)->getQPSolver();
  }

} // namespace CasADi
