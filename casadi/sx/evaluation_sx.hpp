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

#ifndef EVALUATION_SX_HPP
#define EVALUATION_SX_HPP

#include "sx_node.hpp"
#include "output_sx.hpp"
#include "../fx/fx.hpp"
#include <vector>
#include <stack>

namespace CasADi{

/** \brief Represents a function call
  \author Joel Andersson 
  \date 2012
*/
class EvaluationSX : public SXNode{
  private:
    
    /** \brief  Constructor is private, use "create" below (binary version) */
    EvaluationSX(const FX& f, const std::vector<SX>& dep) : f_(f), dep_(dep){}
    
  public:
    
    /** \brief  Create a binary expression */
    inline static std::vector<SX> create(const FX& f, const std::vector<SX>& dep){
      // Assert the right number of inputs
      casadi_assert(f.getNumScalarInputs()==dep.size());
      
      // Allocate outputs
      std::vector<SX> ret(f.getNumScalarOutputs());
      casadi_assert(!ret.empty());
      
      // Create evaluation node
      ret[0] = SX::create(new EvaluationSX(f,dep));
      
      // Create auxillary outputs
      for(int oind=1; oind<ret.size(); ++oind){
	ret[oind] = OutputSX::create(ret[0],oind);
      }
      
      return ret;
    }
    
    /** \brief Destructor
    This might need fixing to avoid stack overflow due to recursive calling.
    */
    virtual ~EvaluationSX(){}
    
    virtual bool isSmooth() const{ return false;}
    
    virtual bool hasDep() const{ return true; }
    
    /** \brief  Number of dependencies */
    virtual int ndep() const{ return dep_.size();}
    
    /** \brief  get the reference of a dependency */
    virtual const SX& dep(int i) const{ return dep_.at(i);}
    virtual SX& dep(int i){ return dep_.at(i);}
    
    /** \brief  Get the operation */
    virtual int getOp() const{ return OP_CALL;}
    
    /** \brief  Print the expression (recursively with a maximum number of levels) */
    virtual void print(std::ostream &stream, long& remaining_calls) const{
      
      // Print the prefix
      stream << f_ << ".call(";
      
      // Print the dependencies
      for(std::vector<SX>::const_iterator i=dep_.begin(); i!=dep_.end(); ++i){
	// Print separator
	if(i!=dep_.begin()) stream << ",";
	
	// Print dependency
	i->print(stream,remaining_calls);
      }
      
      // Print the suffix
      stream << ")";
    }
    
    /** \brief  Function to be called */
    FX f_;
    
    /** \brief  The dependencies of the node (all nonzeros stacked) */
    std::vector<SX> dep_;
};

} // namespace CasADi


#endif // EVALUATION_SX_HPP
