/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include <coin/IpIpoptApplication.hpp>
namespace CasADi{

IpoptSolver::IpoptSolver(){
}
  
IpoptSolver::IpoptSolver(const FX& F, const FX& G, const FX& H, const FX& J){
  assignNode(new IpoptInternal(F,G,H,J));
}

IpoptInternal* IpoptSolver::operator->(){
  return (IpoptInternal*)(NLPSolver::operator->());
}

const IpoptInternal* IpoptSolver::operator->() const{
  return (const IpoptInternal*)(NLPSolver::operator->());
}
    
void IpoptSolver::assertNode() const{
  if(!dynamic_cast<const IpoptInternal*>(get()))
    throw CasadiException("IpoptSolver::assertNode");
}

} // namespace CasADi
