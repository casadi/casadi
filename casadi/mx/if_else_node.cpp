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
  casadi_assert(cond.numel()==1 && cond.size()==1);
  setSparsity(if_true->sparsity());
}

IfNode* IfNode::clone() const{
  return new IfNode(*this);
}

void IfNode::print(std::ostream &stream, const std::vector<std::string>& args) const{
  stream << "(" << args.at(0) << "?" <<  args.at(1) << ":" << args.at(2) << ")";
}

void IfNode::evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
  int nfwd = fwdSens.size();
  int nadj = adjSeed.size();
  bool c = fabs(input[0]->data()[0])<tol;
  vector<double> &outputd = output[0]->data();
  if(c){
    for(int i=0; i<size(); ++i){
      // Function
      outputd[i] = input[0]->data()[i];
      
      // Forward seeds
      for(int d=0; d<nfwd; ++d)
        fwdSens[d][0]->data()[i] = fwdSeed[d][0]->data()[i];
      
      // Adjoint seeds
      for(int d=0; d<nadj; ++d)
        adjSens[d][0]->data()[i] += adjSeed[d][0]->data()[i];
    }
  } else {
    // Zero if false
    for(int i=0; i<size(); ++i){
      outputd[i] = 0;
      for(int d=0; d<nfwd; ++d) 
        fwdSens[d][0]->data()[i] = 0;
    }
  }
}

} // namespace CasADi

