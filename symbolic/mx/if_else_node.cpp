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

#include "if_else_node.hpp"
#include "../stl_vector_tools.hpp"
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

double IfNode::tol = 1e-6;
  
IfNode::IfNode(const MX& cond, const MX& if_true){
  setDependencies(cond,if_true);
  casadi_assert(cond.scalar());
  setSparsity(if_true->sparsity());
}

IfNode* IfNode::clone() const{
  return new IfNode(*this);
}

void IfNode::printPart(std::ostream &stream, int part) const{
  if(part==0){
    stream << "(";
  } else if(part==1){
    stream << " ? ";
  } else {
    stream << ":0)";
  }
}

void IfNode::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  bool cond = fabs(input[0]->data()[0])>tol;
  if(cond){
    
    // Assign if true
    copy(input[1]->begin(),input[1]->end(),output[0]->begin());
    
    // Forward sensitivities
    for(int d=0; d<nfwd; ++d)
      copy(fwdSeed[d][1]->begin(),fwdSeed[d][1]->end(),fwdSens[d][0]->begin());
    
    // Adjoint seeds
    for(int d=0; d<nadj; ++d)
      transform(adjSens[d][1]->begin(),adjSens[d][1]->end(),adjSeed[d][0]->begin(),adjSens[d][1]->begin(),plus<double>());
    
  } else {
    // Zero if false
    fill(output[0]->begin(),output[0]->end(),0.);

    // Forward sensitivities also zero
    for(int d=0; d<nfwd; ++d) 
      fill(fwdSens[d][0]->begin(),fwdSens[d][0]->end(),0.);
  }
}

void IfNode::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
  casadi_error("not implented");
}

void IfNode::evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given){
  casadi_error("not implented");
}

void IfNode::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){ 
  casadi_error("not implented");
}



} // namespace CasADi

