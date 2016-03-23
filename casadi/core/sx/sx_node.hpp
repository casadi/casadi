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
#include "sx_elem.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief  Internal node class for SX
      \author Joel Andersson
      \date 2010
  */
  class SXNode {
    friend class SXElem;
    friend class Matrix<SXElem>;

  public:

    /** \brief  constructor */
    SXNode();

    /** \brief  destructor  */
    virtual ~SXNode();

    ///@{
    /** \brief  check properties of a node */
    virtual bool is_constant() const; // check if constant
    virtual bool is_integer() const; // check if integer
    virtual bool is_symbolic() const; // check if symbolic
    virtual bool hasDep() const; // check if binary
    virtual bool is_zero() const; // check if zero
    virtual bool isAlmostZero(double tol) const; // check if almost zero
    virtual bool is_one() const; // check if one
    virtual bool is_minus_one() const; // check if minus one
    virtual bool isNan() const; // check if not a number
    virtual bool isInf() const; // check if infinity
    virtual bool isMinusInf() const; // check if minus infinity
    ///@}

    ///@{
    /** \brief  Get value of a constant node */
    virtual double to_double() const;  // only works for constant nodes
    virtual int to_int() const;  // only works for integer nodes
    ///@}

    virtual const std::string& name() const; // get the name

    /** \brief get the operation */
    virtual int op() const=0;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const SXNode* node, int depth) const;

    /** \brief  Number of dependencies */
    virtual int ndep() const { return 0;}

    /** \brief  get the reference of a child */
    virtual const SXElem& dep(int i) const;

    /** \brief  get the reference of a child */
    virtual SXElem& dep(int i);

    /** \brief  Check if smooth */
    virtual bool is_smooth() const;

    /** \brief  print */
    virtual void print(std::ostream &stream) const;

    /** \brief Find out which nodes can be inlined */
    void can_inline(std::map<const SXNode*, int>& nodeind) const;

    /** \brief Print compact */
    std::string print_compact(std::map<const SXNode*, int>& nodeind,
                             std::vector<std::string>& intermed) const;

    /** \brief  Print expression */
    virtual std::string print(const std::string& arg1, const std::string& arg2) const = 0;

    // Check if marked (i.e. temporary is negative)
    bool marked() const;

    // Mark by flipping the sign of the temporary and decreasing by one
    void mark();

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
