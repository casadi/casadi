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

#include "scpgen_internal.hpp"

using namespace std;

namespace CasADi{

  SCPgen::SCPgen(){
  }
  
#ifndef WITHOUT_PRE_1_9_X
  SCPgen::SCPgen(const Function& F, const Function& G){
    assignNode(new SCPgenInternal(joinFG(F,G)));
  }
#endif

  SCPgen::SCPgen(const Function& nlp){
    assignNode(new SCPgenInternal(nlp));
  }

  SCPgenInternal* SCPgen::operator->(){
    return static_cast<SCPgenInternal*>(NLPSolver::operator->());
  }

  const SCPgenInternal* SCPgen::operator->() const{
    return static_cast<const SCPgenInternal*>(NLPSolver::operator->());
  }
    
  bool SCPgen::checkNode() const{
    return dynamic_cast<const SCPgenInternal*>(get())!=0;
  }

  const QPSolver SCPgen::getQPSolver() const {
    return (*this)->getQPSolver();
  }
    

} // namespace CasADi
