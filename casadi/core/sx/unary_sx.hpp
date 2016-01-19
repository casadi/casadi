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


#ifndef UNARY_SXElem_HPP
#define UNARY_SXElem_HPP

#include "sx_node.hpp"
#include <stack>

/// \cond INTERNAL

namespace casadi {

/** \brief Represents a basic unary operation on an SXElem node
  \author Joel Andersson
  \date 2012
*/
class CASADI_EXPORT UnarySX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below */
    UnarySX(unsigned char op, const SXElem& dep) : op_(op), dep_(dep) {}

  public:

    /** \brief  Create a unary expression */
    inline static SXElem create(unsigned char op, const SXElem& dep) {
      if (dep.is_constant()) {
        // Evaluate constant
        double dep_val(dep);
        double ret_val;
        casadi_math<double>::fun(op, dep_val, dep_val, ret_val);
        return ret_val;
      } else {
        // Expression containing free variables
        return SXElem::create(new UnarySX(op, dep));
      }
    }

    /** \brief Destructor */
    virtual ~UnarySX() {}

    virtual bool is_smooth() const { return operation_checker<SmoothChecker>(op_);}

    virtual bool hasDep() const { return true; }

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const SXNode* node, int depth) const {
      const UnarySX* n = dynamic_cast<const UnarySX*>(node);
      return n && n->op_ == op_ &&  SXElem::is_equal(n->dep_, dep_, depth-1);
    }

    /** \brief  Number of dependencies */
    virtual int ndep() const { return 1;}

    /** \brief  get the reference of a dependency */
    virtual const SXElem& dep(int i) const { return dep_; }
    virtual SXElem& dep(int i) { return dep_; }

    /** \brief  Get the operation */
    virtual int op() const { return op_;}

    /** \brief  Print expression */
    virtual std::string print(const std::string& arg1, const std::string& arg2) const {
      std::stringstream ss;

      // Print the prefix
      casadi_math<double>::printPre(op_, ss);

      // Print the dependency
      ss << arg1;

      // Print the suffix
      casadi_math<double>::printPost(op_, ss);

      return ss.str();
    }

    /** \brief  The binary operation as an 1 byte integer (allows 256 values) */
    unsigned char op_;

    /** \brief  The dependencies of the node */
    SXElem dep_;
};

} // namespace casadi

/// \endcond
#endif // UNARY_SXElem_HPP
