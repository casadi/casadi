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


#include "sx_node.hpp"
#include <limits>
#include <typeinfo>

using namespace std;
namespace casadi {

  SXNode::SXNode() {
    count = 0;
    temp = 0;
  }

  SXNode::~SXNode() {
    // Make sure that this is there are no scalar expressions pointing to it when it is destroyed
    casadi_assert_warning(count==0,
                          "Reference counting failure. "
                          "Possible cause: Circular dependency in user code.");
  }

  double SXNode::getValue() const {
    return numeric_limits<double>::quiet_NaN();
    /*  std::cerr << "getValue() not defined for class " << typeid(*this).name() << std::endl;
        throw "SXNode::getValue()";*/
  }

  int SXNode::getIntValue() const {
    throw CasadiException(string("getIntValue() not defined for class ") + typeid(*this).name());
  }

  bool SXNode::isZero() const {
    return false;
  }

  bool SXNode::isAlmostZero(double tol) const {
    return false;
  }

  bool SXNode::isOne() const {
    return false;
  }

  bool SXNode::isMinusOne() const {
    return false;
  }

  bool SXNode::isNan() const {
    return false;
  }

  bool SXNode::isInf() const {
    return false;
  }

  bool SXNode::isMinusInf() const {
    return false;
  }

  bool SXNode::isConstant() const {
    return false;
  }

  bool SXNode::isInteger() const {
    return false;
  }

  bool SXNode::isSymbolic() const {
    return false;
  }

  bool SXNode::hasDep() const {
    return false;
  }

  bool SXNode::isEqual(const SXNode* node, int depth) const {
    return false;
  }

  const std::string& SXNode::getName() const {
    throw CasadiException("SXNode::getName failed, the node must be symbolic");
  }

  const SXElement& SXNode::dep(int i) const {
    casadi_error("child() not defined for class " << typeid(*this).name());
  }

  SXElement& SXNode::dep(int i) {
    casadi_error("child() not defined for class " << typeid(*this).name());
  }

  bool SXNode::isSmooth() const {
    return true; // nodes are smooth by default
  }

  void SXNode::print(std::ostream &stream) const {
    long remaining_calls = max_num_calls_in_print_;
    print(stream, remaining_calls);
  }

  bool SXNode::marked() const {
    return temp<0;
  }

  void SXNode::mark() {
    temp = -temp-1;
  }

  long SXNode::max_num_calls_in_print_ = 10000;
  int SXNode::eq_depth_ = 1;

} // namespace casadi
