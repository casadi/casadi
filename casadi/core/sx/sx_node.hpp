/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SX_NODE_HPP
#define CASADI_SX_NODE_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

/** \brief  Scalar expression (which also works as a smart pointer class to this class) */
#include "sx_element.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief  Internal node class for SX
      \author Joel Andersson
      \date 2010
  */
  class SXNode {
    friend class SXElement;
    friend class Matrix<SXElement>;

  public:

    /** \brief  constructor */
    SXNode();

    /** \brief  destructor  */
    virtual ~SXNode();

    ///@{
    /** \brief  check properties of a node */
    virtual bool isConstant() const; // check if constant
    virtual bool isInteger() const; // check if integer
    virtual bool isSymbolic() const; // check if symbolic
    virtual bool hasDep() const; // check if binary
    virtual bool isZero() const; // check if zero
    virtual bool isAlmostZero(double tol) const; // check if almost zero
    virtual bool isOne() const; // check if one
    virtual bool isMinusOne() const; // check if minus one
    virtual bool isNan() const; // check if not a number
    virtual bool isInf() const; // check if infinity
    virtual bool isMinusInf() const; // check if minus infinity
    ///@}

    ///@{
    /** \brief  Get value of a constant node */
    virtual double getValue() const;  // only works for constant nodes
    virtual int getIntValue() const;  // only works for integer nodes
    ///@}

    virtual const std::string& getName() const; // get the name

    /** \brief get the operation */
    virtual int getOp() const=0;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool isEqual(const SXNode* node, int depth) const;

    /** \brief  Number of dependencies */
    virtual int ndep() const { return 0;}

    /** \brief  get the reference of a child */
    virtual const SXElement& dep(int i) const;

    /** \brief  get the reference of a child */
    virtual SXElement& dep(int i);

    /** \brief  Initialize the node (currently used only to give a similar interface to MXNode) */
    void init() {}

    /** \brief  Check if smooth */
    virtual bool isSmooth() const;

    /** \brief  print */
    virtual void print(std::ostream &stream) const;

    /** \brief  print */
    virtual void print(std::ostream &stream, long& remaining_calls) const = 0;

    // Check if marked (i.e. temporary is negative)
    bool marked() const;

    // Mark by flipping the sign of the temporary and decreasing by one
    void mark();

    // Maximum number of calls
    static long max_num_calls_in_print_;

    // Depth when checking equalities
    static int eq_depth_;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    int temp;

    // Reference counter -- counts the number of parents of the node
    unsigned int count;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_SX_NODE_HPP
