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

#include "shared_object_internal.hpp"
#ifdef WITH_EXTRA_CHECKS
#include "function.hpp"
#endif // WITH_EXTRA_CHECKS
#include <typeinfo>

namespace casadi {

  SharedObject::SharedObject() {
    node = nullptr;
  }

  SharedObject::SharedObject(const SharedObject& ref) {
    node = ref.node;
    count_up();
  }

  SharedObject::~SharedObject() {
    count_down();
  }

  void SharedObject::own(SharedObjectInternal* node_) {
    count_down();
    node = node_;
    count_up();
  }

  void SharedObject::assign(SharedObjectInternal* node_) {
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

  SharedObjectInternal* SharedObject::get() const {
    return node;
  }

  bool SharedObject::is_null() const {
    return node==nullptr;
  }

  void SharedObject::count_up() {
#ifdef WITH_EXTRA_CHECKS
    casadi_assert_dev(Function::call_depth_==0);
#endif // WITH_EXTRA_CHECKS
    if (node) node->count++;
  }

  void SharedObject::count_down() {
#ifdef WITH_EXTRA_CHECKS
    casadi_assert_dev(Function::call_depth_==0);
#endif // WITH_EXTRA_CHECKS
    if (!node) return;
    if (node->weak_ref_) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      auto mutex = node->weak_ref_->get_mutex();
      // Avoid triggering a delete while a weak_ref.shared_if_alive is being called
      std::lock_guard<std::mutex> lock(*mutex);
      // Could it be that this mutex is destroyed when the lock goes out of scope?
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

      if (--node->count == 0) {
        delete node;
        node = nullptr;
      }
    } else {
      if (--node->count == 0) {
        delete node;
        node = nullptr;
      }
    }
  }

  SharedObjectInternal* SharedObject::operator->() const {
    casadi_assert_dev(!is_null());
    return node;
  }

  std::string SharedObject::class_name() const {
    return (*this)->class_name();
  }

  void SharedObject::disp(std::ostream& stream, bool more) const {
    if (is_null()) {
      stream << "NULL";
    } else {
      (*this)->disp(stream, more);
    }
  }

  void SharedObject::print_ptr(std::ostream &stream) const {
    stream << node;
  }

  void SharedObject::swap(SharedObject& other) {
    SharedObject temp = *this;
    *this = other;
    other = temp;
  }

  casadi_int SharedObject::getCount() const {
    return (*this)->getCount();
  }

  WeakRef* SharedObject::weak() {
    return (*this)->weak();
  }

  casadi_int SharedObject::__hash__() const {
    return reinterpret_cast<casadi_int>(get());
  }

  WeakRef::WeakRef(int dummy) {
    casadi_assert_dev(dummy==0);
  }

  bool WeakRef::alive() const {
    return !is_null() && (*this)->raw_ != nullptr;
  }

  SharedObject WeakRef::shared() const {
    SharedObject ret;
    if (alive()) {
      ret.own((*this)->raw_);
    }
    return ret;
  }

  bool WeakRef::shared_if_alive(SharedObject& shared) const {
    if (is_null()) return false;
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    // Safe access to ...
    std::lock_guard<std::mutex> lock(*(*this)->mutex_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
    if (alive()) {
      shared.own((*this)->raw_);
      return true;
    }
    return false;
  }

  const WeakRefInternal* WeakRef::operator->() const {
    return static_cast<const WeakRefInternal*>(SharedObject::operator->());
  }

  WeakRefInternal* WeakRef::operator->() {
    return static_cast<WeakRefInternal*>(SharedObject::operator->());
  }

  WeakRef::WeakRef(SharedObject shared) {
    own(shared.weak()->get());
  }

  WeakRef::WeakRef(SharedObjectInternal* raw) {
    own(new WeakRefInternal(raw));
  }

  void WeakRef::kill() {
    casadi_assert_dev((*this)->raw_);
    (*this)->raw_ = nullptr;
  }

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::shared_ptr<std::mutex> WeakRef::get_mutex() const {
    return (*this)->mutex_;
  }
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

} // namespace casadi
