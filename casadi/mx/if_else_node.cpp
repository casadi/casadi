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
#include <cassert>
#include <vector>
#include <sstream>

using namespace std;

namespace CasADi{

double IfElseNode::tol = 1e-6;
  
IfElseNode::IfElseNode(const MX& cond, const MX& if_true, const MX& if_false){
  setDependencies(cond,if_true,if_false);
  assert(cond.numel()==1);
  assert(if_true.size1()==if_false.size1() || if_true.size2()==if_false.size2());
  setSize(if_true.size1(),if_true.size2());
}

IfElseNode* IfElseNode::clone() const{
  return new IfElseNode(*this);
}

void IfElseNode::print(std::ostream &stream) const{
  stream << "(" << dep(0) << "?" <<  dep(1) << ":" << dep(2) << ")";
}

void IfElseNode::evaluate(int fsens_order, int asens_order){
  const vector<double>& cond = input(0);  // condition
  const vector<double>& if_true = input(1);  // if condition true
  const vector<double>& if_false = input(2);  // if condition false
  vector<double>& res = output();
  if(fabs(cond[0])>tol)
    copy(if_true.begin(),if_true.end(),res.begin());
  else
    copy(if_false.begin(),if_false.end(),res.begin());
  
  if(fsens_order>0){
    const vector<double>& if_true_der = fwdSeed(1);  // if condition true derivative
    const vector<double>& if_false_der = fwdSeed(2);  // if condition false derivative
    vector<double>& res_der = fwdSens();
    if(fabs(cond[0])>tol)
      copy(if_true_der.begin(),if_true_der.end(),res_der.begin());
    else
      copy(if_false_der.begin(),if_false_der.end(),res_der.begin());
  }

  if(asens_order>0){
    vector<double>& if_true_der = adjSens(1);  // if condition true derivative
    vector<double>& if_false_der = adjSens(2);  // if condition false derivative
    const vector<double>& res_der = adjSeed();
    if(fabs(cond[0])>tol){ // if true
      for(int i=0; i<res_der.size(); ++i)
        if_true_der[i] += res_der[i];
    } else { // if false
      for(int i=0; i<res_der.size(); ++i)
        if_false_der[i] += res_der[i];
    }
  }
}

} // namespace CasADi

