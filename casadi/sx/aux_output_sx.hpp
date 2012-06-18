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

#ifndef AUX_OUTPUT_SX_HPP
#define AUX_OUTPUT_SX_HPP

#include "sx_node.hpp"

namespace CasADi{

/** \brief Represents an auxillary output of a function
  \author Joel Andersson 
  \date 2012
*/
class AuxOutputSX : public SXNode{
  private:
    
    /** \brief  Constructor is private, use "create" below (binary version) */
    AuxOutputSX(const SX& dep, int oind) : dep_(dep), oind_(oind){}
    
  public:
    
    /** \brief  Create an auxillary output */
    inline static SX create(const SX& dep, int oind){
      return SX::create(new AuxOutputSX(dep,oind));
    }
    
    /** \brief Destructor
    This might need fixing to avoid stack overflow due to recursive calling.
    */
    virtual ~AuxOutputSX(){}
    
    virtual bool isSmooth() const{ return false;}
    
    virtual bool hasDep() const{ return true; }
    
    /** \brief  Number of dependencies */
    virtual int ndep() const{ return 1;}
    
    /** \brief  get the reference of a dependency */
    virtual const SX& dep(int i) const{ return dep_;}
    virtual SX& dep(int i){ return dep_;}
    
    /** \brief  Get the operation */
    virtual int getOp() const{ return OP_FOUTPUT;}
    
    /** \brief  Print the expression (recursively with a maximum number of levels) */
    virtual void print(std::ostream &stream, long& remaining_calls) const{

      // Print dependency
      dep_->print(stream,remaining_calls);
      
      // Print the index
      stream << "[" << oind_ << "]";
    }
      
    
    /** \brief The main function which includes the first output */
    SX dep_;
    
    /** \brief Output index */
    int oind_;
};

 
} // namespace CasADi


#endif // AUX_OUTPUT_SX_HPP
