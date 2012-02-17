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

#ifndef BINARY_SCALAR_HPP
#define BINARY_SCALAR_HPP

#include "sx_node.hpp"
#include <stack>

namespace CasADi{

/** \brief Represents a basic binary operation on two SX nodes
  \author Joel Andersson 
  \date 2010
*/
class BinarySXNode : public SXNode{
  private:
    
    /** \brief  Constructor is private, use "create" below (unary version) */
    BinarySXNode(unsigned char op_, const SX& child1_){
      op = op_;
      child[0] = child1_;
      child[1] = casadi_limits<SX>::zero;
    }
    
    /** \brief  Constructor is private, use "create" below (binary version) */
    BinarySXNode(unsigned char op_, const SX& child1_, const SX& child2_){
      op = op_;
      child[0] = child1_;
      child[1] = child2_;
    }
    
  public:
    
    /** \brief  Create a unary expression */
    inline static SX create(unsigned char op_, const SX& child1_){
      return SX::create(new BinarySXNode(op_,child1_));
    }
    
    /** \brief  Create a binary expression */
    inline static SX create(unsigned char op_, const SX& child1_, const SX& child2_){
      return SX::create(new BinarySXNode(op_,child1_,child2_));
    }
    
    /** \brief  Create a unary expression (templated version) */
    template<int op_>
    inline static SX createT(const SX& child1_){
      return SX::create(new BinarySXNode(op_,child1_));
    }
    
    /** \brief  Create a binary expression (templated version) */
    template<int op_>
    inline static SX createT(const SX& child1_, const SX& child2_){
      return SX::create(new BinarySXNode(op_,child1_,child2_));
    }
    
    /** \brief Destructor
    This is a rather complex destructor which is necessary since the default destructor 
    can cause stack overflow due to recursive calling.
    */
    virtual ~BinarySXNode(){
      // Start destruction method if any of the dependencies has dependencies
      for(int c1=0; c1<2; ++c1){
        // Get the node of the child and remove it from the smart pointer
        SXNode* n1 = child[c1].assignNoDelete(casadi_limits<SX>::nan);
        
        // Check if this was the last reference
        if(n1->count==0){

          // Check if binary
          if(!n1->hasDep()){ // n1 is not binary

            delete n1; // Delete stright away 

          } else { // n1 is binary
            
            // Stack of experssions to be deleted
            std::stack<SXNode*> deletion_stack;
            
            // Add the node to the deletion stack
            deletion_stack.push(n1);
            
            // Process stack
            while(!deletion_stack.empty()){
              
              // Top element
              SXNode *t = deletion_stack.top();
              
              // Check if the top element has dependencies with dependencies
              bool added_to_stack = false;
              for(int c2=0; c2<2; ++c2){ // for all children of the child
                
                // Get the node of the child of the top element and remove it from the smart pointer
                SXNode *n2 = t->dep(c2).assignNoDelete(casadi_limits<SX>::nan);
                
                // Check if this is the only reference to the element
                if(n2->count == 0){
                  
                  // Check if binary
                  if(!n2->hasDep()){
                    
                    // Delete stright away if not binary
                    delete n2;
                    
                  } else {
                    
                    // Add to deletion stack
                    deletion_stack.push(n2);
                    added_to_stack = true;
                  }
                }
              }
              
              // Delete and pop from stack if nothing added to the stack
              if(!added_to_stack){
                delete deletion_stack.top();
                deletion_stack.pop();
              }
            } // while
          }
        }
      }
    }
    
    virtual bool isSmooth() const{ return operation_checker<SmoothChecker>(op);}
    
    virtual bool hasDep() const{ return true; }
    
    /** \brief  Number of dependencies */
    virtual int ndep() const{ return 2;}
    
    /** \brief  get the reference of a child */
    virtual const SX& dep(int i) const{ return child[i];}
    virtual SX& dep(int i){ return child[i];}
    
    /** \brief  Get the operation */
    virtual int getOp() const{ return op;}
    
    /** \brief  Print the expression (recursively with a maximum number of levels) */
    virtual void print(std::ostream &stream, long& remaining_calls) const{

      // Print the prefix
      casadi_math<double>::printPre(op,stream);
      
      // Print the first child
      child[0].print(stream,remaining_calls);
      
      // If binary
      if(casadi_math<double>::ndeps(op)>1){
        
        //Print the infix
        casadi_math<double>::printSep(op,stream);

        // Print the second child
        child[1].print(stream,remaining_calls);
      }
      
      // Print the suffix
      casadi_math<double>::printPost(op,stream);
    }
    
    /** \brief  The binary operation as an 1 byte integer (allows 256 values) */
    unsigned char op;
    
    /** \brief  The children of the node */
    SX child[2];
};

} // namespace CasADi


#endif // BINARY_SCALAR_HPP
