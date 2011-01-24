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

#include "binary_sx_node.hpp"
#include <cassert>
#include <stack>

using namespace std;
namespace CasADi{

BinarySXNode::~BinarySXNode(){
  for(int c=0; c<2; ++c){
    if(child[c]->isBinary()){
      if(--child[c].node->count == 0){
        // Do not delete directly, add to a stack instead
        stack<BinarySXNode*> s;
        s.push((BinarySXNode*)child[c].node);

        // Process the stack
        while(!s.empty()){
          BinarySXNode* t = s.top();
          // Delete the nodes of both children
          for(int cc=0; cc<2; ++cc){
            if(t->child[cc]->isBinary()){
              if(--t->child[cc].node->count == 0){ // "delete"
                // Add to delete stack
                s.push((BinarySXNode*)t->child[cc].node);
              }
              // Replace the binary node with a nan-node
              t->child[cc].node = casadi_limits<SX>::nan.node;
              t->child[cc].node->count++;
            }
          }
          
          // If the node is still on the top, delete it
          if(t==s.top()){
            delete t;
            s.pop();
          }
        }
      }
      
      // Replace the binary node with a nan-node
      child[c].node = casadi_limits<SX>::nan.node;
      child[c].node->count++;
    }
  }
}

  
void BinarySXNode::print(ostream &stream) const{
  stringstream s0,s1;
  s0 << child[0];
  s1 << child[1];
  print_c[op](stream,s0.str(),s1.str());
}

bool BinarySXNode::isSmooth() const{
  if(op == STEP_NODE || op == FLOOR_NODE)
    return false;
  else
    return true;
}

const SX& BinarySXNode::dependent(int i) const{
  assert(i==0 || i==1);
  return child[i];
}

} // namespace CasADi
