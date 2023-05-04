/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include <math.h>
#include <sstream>
#include <string>

/** \brief  Scalar expression (which also works as a smart pointer class to this class)

    \identifier{9s} */
#include "sx_elem.hpp"


/// \cond INTERNAL
namespace casadi {

  /** \brief  Internal node class for SX

      \author Joel Andersson
      \date 2010

      \identifier{9t} */
  class SXNode {
    friend class SXElem;
    friend class Matrix<SXElem>;
    friend class UniversalNodeOwner;
  public:

    /** \brief  constructor

        \identifier{9u} */
    SXNode();

    /** \brief  destructor

        \identifier{9v} */
    virtual ~SXNode();

    ///@{
    /** \brief  check properties of a node

        \identifier{9w} */
    virtual bool is_constant() const { return false; }
    virtual bool is_integer() const { return false; }
    virtual bool is_symbolic() const { return false; }
    virtual bool is_zero() const { return false; }
    virtual bool is_op(casadi_int op) const { return false; }
    virtual bool is_almost_zero(double tol) const { return false; }
    virtual bool is_one() const { return false; }
    virtual bool is_minus_one() const { return false; }
    virtual bool is_nan() const { return false; }
    virtual bool is_inf() const { return false; }
    virtual bool is_minus_inf() const { return false; }
    ///@}

    ///@{
    /** \brief  Get value of a constant node

        \identifier{9x} */
    virtual double to_double() const;  // only works for constant nodes
    virtual casadi_int to_int() const;  // only works for integer nodes
    ///@}

    // get the name
    virtual const std::string& name() const;

    /** \brief Get type name

        \identifier{9y} */
    virtual std::string class_name() const = 0;

    /** \brief get the operation

        \identifier{9z} */
    virtual casadi_int op() const=0;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{a0} */
    virtual bool is_equal(const SXNode* node, casadi_int depth) const;

    /** \brief  Number of dependencies

        \identifier{a1} */
    virtual casadi_int n_dep() const { return 0;}

    /** \brief  get the reference of a child

        \identifier{a2} */
    virtual const SXElem& dep(casadi_int i) const;

    /** \brief  get the reference of a child

        \identifier{a3} */
    virtual SXElem& dep(casadi_int i);

    /** \brief  Check if smooth

        \identifier{a4} */
    virtual bool is_smooth() const { return true; }

    /** \brief  print

        \identifier{a5} */
    virtual void disp(std::ostream& stream, bool more) const;

    /** \brief Find out which nodes can be inlined

        \identifier{a6} */
    void can_inline(std::map<const SXNode*, casadi_int>& nodeind) const;

    /** \brief Print compact

        \identifier{a7} */
    std::string print_compact(std::map<const SXNode*, casadi_int>& nodeind,
                             std::vector<std::string>& intermed) const;

    /** \brief  Print expression

        \identifier{a8} */
    virtual std::string print(const std::string& arg1, const std::string& arg2) const = 0;

    // Check if marked (i.e. temporary is negative)
    bool marked() const;

    // Mark by flipping the sign of the temporary and decreasing by one
    void mark() const;

    /** \brief Non-recursive delete

        \identifier{a9} */
    static void safe_delete(SXNode* n);

    // Depth when checking equalities
    static casadi_int eq_depth_;

    /** Temporary variables to be used in user algorithms like sorting,
        the user is responsible of making sure that use is thread-safe
        The variable is initialized to zero
    */
    mutable int temp;

    // Reference counter -- counts the number of parents of the node
    unsigned int count;

    /** \brief Serialize an object

        \identifier{aa} */
    void serialize(SerializingStream& s) const;

    virtual void serialize_node(SerializingStream& s) const;

    static SXNode* deserialize(DeserializingStream& s);

    static std::map<casadi_int, SXNode* (*)(DeserializingStream&)> deserialize_map;


  };

} // namespace casadi
/// \endcond
#endif // CASADI_SX_NODE_HPP
