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
  // Start destruction method if any of the dependencies has dependencies
  for(int c=0; c<2; ++c){
    if(child[c]->hasDep()){
      // Check if there are other "owners" of the node
      if(child[c]->count != 1){
        
        // Replace with a nan-node
        child[c] = casadi_limits<SX>::nan;
      } else {
      
        // Create stack of nodes to be deleted
        stack<SX> s;

        // Add the child
        s.push(child[c]);
        child[c] = casadi_limits<SX>::nan;
        
        // Process stack
        while(!s.empty()){
          // Top element
          SX t = s.top();
          
          // Check if the top element has dependencies with dependencies
          bool found_dep = false;
          for(int i=0; i<t->ndep(); ++i){
            if(t->dep(i)->hasDep()){
              // Check if this is the only reference to the element
              if(t->dep(i)->count==1){
                // Remove and add to stack
                s.push(t->dep(i));
                t->dep(i) = casadi_limits<SX>::nan;
                found_dep = true;
                break;
              } else {
                // Replace with an element without dependencies
                t->dep(i) = casadi_limits<SX>::nan;
              }
            }
          }
          
          // Pop from stack if no dependencies found
          if(!found_dep){
            s.pop();
          }
        }
      }
    }
  }
}
  
void BinarySXNode::print(ostream &stream) const{
  casadi_math<double>::printPre[op](stream);
  stream << child[0];
  if (casadi_math<double>::ndeps[op]>1) {
    casadi_math<double>::printSep[op](stream);
    stream << child[1];
  }
  casadi_math<double>::printPost[op](stream);
}

bool BinarySXNode::isSmooth() const{
  if(op == STEP || op == FLOOR)
    return false;
  else
    return true;
}

const SX& BinarySXNode::dep(int i) const{
  assert(i==0 || i==1);
  return child[i];
}

SX& BinarySXNode::dep(int i){
  assert(i==0 || i==1);
  return child[i];
}

} // namespace CasADi
