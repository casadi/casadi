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

#include "../../symbolic/stl_vector_tools.hpp"
#include "../../symbolic/matrix/matrix_tools.hpp"

using namespace std;
namespace CasADi {

DSDPInternal* DSDPInternal::clone() const{
  // Return a deep copy
  DSDPInternal* node = new DSDPInternal(input(SDP_C).sparsity(),input(SDP_A).sparsity());
  if(!node->is_init_)
    node->init();
  return node;
}
  
DSDPInternal::DSDPInternal(const CRSSparsity &C, const CRSSparsity &A) : SDPSolverInternal(C,A){
  
}

DSDPInternal::~DSDPInternal(){ 

}

void DSDPInternal::init(){
  SDPSolverInternal::init();
}

void DSDPInternal::evaluate(int nfdir, int nadir) {
  
}

} // namespace CasADi
