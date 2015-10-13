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


#include "shared_object.hpp"
#include "casadi_exception.hpp"
#include "weak_ref.hpp"

#include <typeinfo>

using namespace std;
namespace casadi {

  SharedObject::SharedObject() {
    node = 0;
  }

  SharedObjectNode::SharedObjectNode(const SharedObjectNode& node) {
    count = 0; // reference counter is _not_ copied
    weak_ref_ = 0; // nor will they have the same weak references
  }

  SharedObjectNode& SharedObjectNode::operator=(const SharedObjectNode& node) {
    // do _not_ copy the reference counter
    return *this;
  }

  SharedObject::SharedObject(const SharedObject& ref) {
    node = ref.node;
    count_up();
  }

  SharedObject::~SharedObject() {
    count_down();
  }

  void SharedObject::assignNode(SharedObjectNode* node_) {
    count_down();
    node = node_;
    count_up();
  }

  void SharedObject::assignNodeNoCount(SharedObjectNode* node_) {
    node = node_;
  }

  SharedObject& SharedObject::operator=(const SharedObject& ref) {
    // quick return if the old and new pointers point to the same object
    if (node == ref.node) return *this;

    // decrease the counter and delete if this was the last pointer
    count_down();

    // save the new pointer
    node = ref.node;
    count_up();
    return *this;
  }

  const SharedObjectNode* SharedObject::get() const {
    return node;
  }

  SharedObjectNode* SharedObject::get() {
    return node;
  }

  bool SharedObject::isNull() const {
    return node==0;
  }

  void SharedObject::count_up() {
    if (node) node->count++;
  }

  void SharedObject::count_down() {
    if (node && --node->count == 0) {
      delete node;
      node = 0;
    }
  }

  const SharedObjectNode* SharedObject::operator->() const {
    casadi_assert(!isNull());
    return node;
  }

  SharedObjectNode* SharedObject::operator->() {
    casadi_assert(!isNull());
    return node;
  }

  SharedObjectNode::SharedObjectNode() {
    count = 0;
    weak_ref_ = 0;
  }

  SharedObjectNode::~SharedObjectNode() {
    if (count!=0) {
      // Note that casadi_assert_warning cannot be used in destructors
      std::cerr << "Reference counting failure." <<
                   "Possible cause: Circular dependency in user code." << std::endl;
    }
    if (weak_ref_!=0) {
      weak_ref_->kill();
      delete weak_ref_;
    }
  }

  void SharedObject::repr(std::ostream &stream, bool trailing_newline) const {
    if (isNull()) {
      stream << 0;
    } else {
      (*this)->repr(stream);
    }
    if (trailing_newline) stream << std::endl;
  }

  void SharedObjectNode::repr(std::ostream &stream) const {
    // Print description by default
    print(stream);
  }

  void SharedObject::print(std::ostream &stream, bool trailing_newline) const {
    if (isNull()) {
      stream << "Null pointer of class \"" << typeid(this).name() << "\"";
    } else {
      (*this)->print(stream);
    }
    if (trailing_newline) stream << std::endl;
  }

  void SharedObject::printPtr(std::ostream &stream) const {
    stream << node;
  }

  void SharedObjectNode::print(std::ostream &stream) const {
    // Print the name of the object by default
    stream << typeid(this).name();
  }

  void SharedObject::swap(SharedObject& other) {
    SharedObject temp = *this;
    *this = other;
    other = temp;
  }

  int SharedObject::getCount() const {
    return (*this)->getCount();
  }

  int SharedObjectNode::getCount() const {
    return count;
  }

  WeakRef* SharedObject::weak() {
    return (*this)->weak();
  }

  WeakRef* SharedObjectNode::weak() {
    if (weak_ref_==0) {
      weak_ref_ = new WeakRef(this);
    }
    return weak_ref_;
  }

  size_t SharedObject::__hash__() const {
    return reinterpret_cast<size_t>(get());
  }

} // namespace casadi


