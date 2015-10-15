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
    if (count!=0) {
      // Note that casadi_assert_warning cannot be used in destructors
      std::cerr << "Reference counting failure." <<
                   "Possible cause: Circular dependency in user code." << std::endl;
    }
  }

  double SXNode::getValue() const {
    return numeric_limits<double>::quiet_NaN();
    /*  userOut<true, PL_WARN>() << "getValue() not defined for class " << typeid(*this).name() << std::endl;
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

  bool SXNode::zz_is_equal(const SXNode* node, int depth) const {
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
    // Find out which noded can be inlined
    std::map<const SXNode*, int> nodeind;
    can_inline(nodeind);

    // Print expression
    vector<string> intermed;
    string s = printCompact(nodeind, intermed);

    // Print intermediate expressions
    for (int i=0; i<intermed.size(); ++i)
      stream << "@" << (i+1) << "=" << intermed[i] << ", ";

    // Print this
    stream << s;
  }

  bool SXNode::marked() const {
    return temp<0;
  }

  void SXNode::mark() {
    temp = -temp-1;
  }

  void SXNode::can_inline(std::map<const SXNode*, int>& nodeind) const {
    // Add or mark node in map
    std::map<const SXNode*, int>::iterator it=nodeind.find(this);
    if (it==nodeind.end()) {
      // First time encountered, mark inlined
      nodeind.insert(it, make_pair(this, 0));

      // Handle dependencies with recursion
      for (int i=0; i<ndep(); ++i) {
        dep(i)->can_inline(nodeind);
      }
    } else if (it->second==0 && getOp()!=OP_PARAMETER) {
      // Node encountered before, do not inline (except if symbolic primitive)
      it->second = -1;
    }
  }

  std::string SXNode::printCompact(std::map<const SXNode*, int>& nodeind,
                                   std::vector<std::string>& intermed) const {
    // Get reference to node index
    int& ind = nodeind[this];

    // If positive, already in intermediate expressions
    if (ind>0) {
      stringstream ss;
      ss << "@" << ind;
      return ss.str();
    }

    // Get expressions for dependencies
    std::string arg[2];
    for (int i=0; i<ndep(); ++i) {
      arg[i] = dep(i)->printCompact(nodeind, intermed);
    }

    // Get expression for this
    string s = print(arg[0], arg[1]);

    // Decide what to do with the expression
    if (ind==0) {
      // Inline expression
      return s;
    } else {
      // Add to list of intermediate expressions and return reference
      intermed.push_back(s);
      ind = intermed.size(); // For subsequent references
      stringstream ss;
      ss << "@" << ind;
      return ss.str();
    }
  }

  int SXNode::eq_depth_ = 1;

} // namespace casadi
