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

#include "ipopt_internal.hpp"
#include "ipopt_nlp.hpp"

using namespace std;
namespace CasADi{

  IpoptSolver::IpoptSolver(){
  }
  
  IpoptSolver::IpoptSolver(const FX& F, const FX& G){
    assignNode(new IpoptInternal(joinFG(F,G)));
  }

  IpoptSolver::IpoptSolver(const FX& nlp){
    assignNode(new IpoptInternal(nlp));
  }

  IpoptInternal* IpoptSolver::operator->(){
    return static_cast<IpoptInternal*>(NLPSolver::operator->());
  }

  const IpoptInternal* IpoptSolver::operator->() const{
    return static_cast<const IpoptInternal*>(NLPSolver::operator->());
  }
    
  bool IpoptSolver::checkNode() const{
    return dynamic_cast<const IpoptInternal*>(get());
  }

  DMatrix IpoptSolver::getReducedHessian(){
    return (*this)->getReducedHessian();
  }

} // namespace CasADi
