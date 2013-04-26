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

#ifndef UNARY_SX_HPP
#define UNARY_SX_HPP

#include "sx_node.hpp"
#include <stack>

namespace CasADi{

/** \brief Represents a basic unary operation on an SX node
  \author Joel Andersson 
  \date 2012
*/
class UnarySX : public SXNode{
  private:
    
    /** \brief  Constructor is private, use "create" below */
    UnarySX(unsigned char op, const SX& dep) : op_(op), dep_(dep){}
    
  public:
    
    /** \brief  Create a unary expression */
    inline static SX create(unsigned char op, const SX& dep){
      if(dep.isConstant()){
        // Evaluate constant
        double dep_val = dep.getValue();
        double ret_val;
        casadi_math<double>::fun(op,dep_val,dep_val,ret_val);
        return ret_val;
      } else {
        // Expression containing free variables
        return SX::create(new UnarySX(op,dep));
      }
    }
    
    /** \brief Destructor */
    virtual ~UnarySX(){}
    
    virtual bool isSmooth() const{ return operation_checker<SmoothChecker>(op_);}
    
    virtual bool hasDep() const{ return true; }
    
    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const SXNode* node, int depth) const{
      const UnarySX* n = dynamic_cast<const UnarySX*>(node);
      return n && n->op_ == op_ &&  n->dep_.isEqual(dep_,depth-1);
    }
      
    /** \brief  Number of dependencies */
    virtual int ndep() const{ return 1;}
    
    /** \brief  get the reference of a dependency */
    virtual const SX& dep(int i) const{ return dep_; }
    virtual SX& dep(int i){ return dep_; }
    
    /** \brief  Get the operation */
    virtual int getOp() const{ return op_;}
    
    /** \brief  Print the expression (recursively with a maximum number of levels) */
    virtual void print(std::ostream &stream, long& remaining_calls) const{

      // Print the prefix
      casadi_math<double>::printPre(op_,stream);
      
      // Print the dependency
      dep_.print(stream,remaining_calls);
      
      // Print the suffix
      casadi_math<double>::printPost(op_,stream);
    }
    
    /** \brief  The binary operation as an 1 byte integer (allows 256 values) */
    unsigned char op_;
    
    /** \brief  The dependencies of the node */
    SX dep_;
};

} // namespace CasADi


#endif // UNARY_SX_HPP
