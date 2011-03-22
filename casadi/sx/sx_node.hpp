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

#ifndef SX_NODE_HPP
#define SX_NODE_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

/** \brief  Scalar expression (which also works as a smart pointer class to this class) */
#include "sx.hpp"

namespace CasADi{

/** \brief  Interal node class for SX
  \author Joel Andersson 
  \date 2010
*/
class SXNode{
friend class SX;

public:

/** \brief  constructor */
SXNode();

/** \brief  destructor  */
virtual ~SXNode();
//@{
/** \brief  check properties of a node */
virtual bool isConstant() const; // check if constant
virtual bool isInteger() const; // check if integer
virtual bool isSymbolic() const; // check if symbolic
virtual bool hasDep() const; // check if binary
virtual bool isZero() const; // check if zero
virtual bool isOne() const; // check if one
virtual bool isMinusOne() const; // check if minus one
virtual bool isNan() const; // check if not a number
virtual bool isInf() const; // check if infinity
virtual bool isMinusInf() const; // check if minus infinity
//@}

//@{
/** \brief  Get value of a constant node */
virtual double getValue() const;  // only works for constant nodes
virtual int getIntValue() const;  // only works for integer nodes
//@}

virtual const std::string& getName() const; // get the name
/** \brief get the operation 
only for binary nodes
*/
virtual int getOp() const; // get the operation (only for binary nodes)
/// comparison
bool isEqual(const SXNode& node) const; // comparison
/// comparison
bool isEqual(const SX& scalar) const; // comparison

/** \brief  Number of dependencies */
virtual int ndep() const{ return 0;}

/** \brief  get the reference of a child */
virtual const SX& dep(int i) const;

/** \brief  get the reference of a child */
virtual SX& dep(int i);

/** \brief  Check if smooth */
virtual bool isSmooth() const;

/** \brief  print */
virtual void print(std::ostream &stream) const = 0;

/** Temporary variables to be used in user algorithms like sorting, 
 the user is resposible of making sure that use is thread-safe
 The variable is initialized to zero
*/
int temp;
// int temp2;

// Reference counter -- counts the number of parents of the node
unsigned int count;

};

} // namespace CasADi


/** \brief  Derived classes */
#include "constant_sx_node.hpp"
#include "symbolic_sx_node.hpp"
#include "binary_sx_node.hpp"

#endif // SX_NODE_HPP
