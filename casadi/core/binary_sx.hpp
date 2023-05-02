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


#ifndef CASADI_BINARY_SX_HPP
#define CASADI_BINARY_SX_HPP

#include "sx_node.hpp"
#include "serializing_stream.hpp"

/// \cond INTERNAL
namespace casadi {

/** \brief Represents a basic binary operation on two SXElem nodes

  \author Joel Andersson
  \date 2010

    \identifier{115} */
class BinarySX : public SXNode {
  private:

    /** \brief  Constructor is private, use "create" below

        \identifier{116} */
    BinarySX(unsigned char op, const SXElem& dep0, const SXElem& dep1) :
        op_(op), dep0_(dep0), dep1_(dep1) {}

  public:

    /** \brief  Create a binary expression

        \identifier{117} */
    inline static SXElem create(unsigned char op, const SXElem& dep0, const SXElem& dep1) {
      if (dep0.is_constant() && dep1.is_constant()) {
        // Evaluate constant
        double dep0_val(dep0);
        double dep1_val(dep1);
        double ret_val;
        casadi_math<double>::fun(op, dep0_val, dep1_val, ret_val);
        return ret_val;
      } else {
        // Expression containing free variables
        return SXElem::create(new BinarySX(op, dep0, dep1));
      }
    }

    /** \brief Destructor

    This is a rather complex destructor which is necessary since the default destructor
    can cause stack overflow due to recursive calling.

        \identifier{118} */
    ~BinarySX() override {
      safe_delete(dep0_.assignNoDelete(casadi_limits<SXElem>::nan));
      safe_delete(dep1_.assignNoDelete(casadi_limits<SXElem>::nan));
    }

    // Class name
    std::string class_name() const override {return "BinarySX";}

    bool is_smooth() const override { return operation_checker<SmoothChecker>(op_);}

    bool is_op(casadi_int op) const override { return op_==op; }

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{119} */
    bool is_equal(const SXNode* node, casadi_int depth) const override {
      const BinarySX* n = dynamic_cast<const BinarySX*>(node);
      if (n==nullptr) return false;
      if (n->op_ != op_) return false;
      if (SXElem::is_equal(n->dep0_, dep0_, depth-1)
          && SXElem::is_equal(n->dep1_, dep1_, depth-1)) return true;
      if (operation_checker<CommChecker>(op_)
          && SXElem::is_equal(n->dep1_, dep0_, depth-1)
          && SXElem::is_equal(n->dep0_, dep1_, depth-1)) return true;
      return false;
    }

    /** \brief  Number of dependencies

        \identifier{11a} */
    casadi_int n_dep() const override { return 2;}

    /** \brief  get the reference of a dependency

        \identifier{11b} */
    const SXElem& dep(casadi_int i) const override { return i==0 ? dep0_ : dep1_;}
    SXElem& dep(casadi_int i) override { return i==0 ? dep0_ : dep1_;}

    /** \brief  Get the operation

        \identifier{11c} */
    casadi_int op() const override { return op_;}

    /** \brief  Print expression

        \identifier{11d} */
    std::string print(const std::string& arg1, const std::string& arg2) const override {
      return casadi_math<double>::print(op_, arg1, arg2);
    }

    /** \brief  The binary operation as an 1 byte integer (allows 256 values)

        \identifier{11e} */
    unsigned char op_;

    /** \brief  The dependencies of the node

        \identifier{11f} */
    SXElem dep0_, dep1_;

    void serialize_node(SerializingStream& s) const override {
      s.pack("UnarySX::dep0", dep0_);
      s.pack("UnarySX::dep1", dep1_);
    }

    /** \brief Deserialize without type information

        \identifier{11g} */
    static SXNode* deserialize(DeserializingStream& s, casadi_int op) {
      SXElem dep0, dep1;
      s.unpack("UnarySX::dep0", dep0);
      s.unpack("UnarySX::dep1", dep1);
      return new BinarySX(op, dep0, dep1);
    }
};

} // namespace casadi
/// \endcond

#endif // CASADI_BINARY_SX_HPP
