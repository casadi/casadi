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

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

using namespace std;
namespace CasADi {

NLPImplicitInternal* NLPImplicitInternal::clone() const{
  // Return a deep copy
  NLPImplicitInternal* node = new NLPImplicitInternal(f_,nrhs_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
NLPImplicitInternal::NLPImplicitInternal(const FX& f, int nrhs) : ImplicitFunctionInternal(f,nrhs) {

  addOption("nlp_solver",       OT_NLPSOLVER, GenericType(), "The NLPSolver used to solve the implicit system.");
  addOption("nlp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLPSolver");
  
}

NLPImplicitInternal::~NLPImplicitInternal(){ 
}

void NLPImplicitInternal::evaluate(int nfdir, int nadir) {
  throw CasadiException("NLPImplicitInternal::evaluate() not implemented yet");
}

void NLPImplicitInternal::init(){

  
  ImplicitFunctionInternal::init();

 
}

} // namespace CasADi

