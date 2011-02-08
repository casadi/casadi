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
#include "binary_functions.hpp"

namespace CasADi{

/** \brief Represents a basic binary operation on two SX nodes
  \author Joel Andersson 
  \date 2010
*/
class BinarySXNode : public SXNode{
public:

/** \brief  Constructors */
BinarySXNode(OPERATION op_, const SX& child1_){
 op = op_;
 child[0] = child1_;
 child[1] = 0;
}

BinarySXNode(OPERATION op_, const SX& child1_, const SX& child2_){
 op = op_;
 child[0] = child1_;
 child[1] = child2_;
}

/// This is a rather complex destructor which is necessary since the default destructor can cause stack overflow due to recursive calling
virtual ~BinarySXNode();

virtual bool isSmooth() const;

virtual bool isBinary() const{ return true; }

/** \brief  Number of dependencies */
virtual int ndep() const{ return 2;}

/** \brief  get the reference of a child */
virtual const SX& dep(int i) const;

/** \brief  Get the operation */
virtual int getOp() const{ return op;}

/** \brief  Data members */
OPERATION  op;
SX      child[2];

virtual void print(std::ostream &stream) const;


};

} // namespace CasADi


#endif // BINARY_SCALAR_HPP
