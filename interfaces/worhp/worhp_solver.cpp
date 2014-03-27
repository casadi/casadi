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

#include "worhp_solver.hpp"
#include "worhp_internal.hpp"

using namespace std;

namespace CasADi{

  WorhpSolver::WorhpSolver(){
  }
  
  WorhpSolver::WorhpSolver(const Function& F, const Function& G){
    assignNode(new WorhpInternal(joinFG(F,G)));
  }

  WorhpSolver::WorhpSolver(const Function& nlp){
    assignNode(new WorhpInternal(nlp));
  }

  WorhpInternal* WorhpSolver::operator->(){
    return static_cast<WorhpInternal*>(NLPSolver::operator->());
  }

  const WorhpInternal* WorhpSolver::operator->() const{
    return static_cast<const WorhpInternal*>(NLPSolver::operator->());
  }
    
  bool WorhpSolver::checkNode() const{
    return dynamic_cast<const WorhpInternal*>(get());
  }
  
  void WorhpSolver::setOptionsFromFile(const std::string & file) {
    dynamic_cast<WorhpInternal*>(get())->setOptionsFromFile(file);
  }

} // namespace CasADi
