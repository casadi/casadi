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

#ifndef INNER_PROD_HPP
#define INNER_PROD_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

namespace CasADi{
  /** \brief Matrix inner_prod
      \author Joel Andersson
      \date 2013
  */
  class InnerProd : public MXNode{
  public:

    /// Constructor
    InnerProd(const MX& x, const MX& y);

    /// Clone function
    virtual InnerProd* clone() const{ return new InnerProd(*this);}
      
    /// Destructor
    virtual ~InnerProd(){}
    
    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;
            
    /** \brief Get the operation */
    virtual int getOp() const{ return OP_INNER_PROD;}    
  };


} // namespace CasADi

#endif // INNER_PROD_HPP
