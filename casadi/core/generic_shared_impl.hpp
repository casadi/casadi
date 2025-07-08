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


#ifndef CASADI_GENERIC_SHARED_IMPL_HPP
#define CASADI_GENERIC_SHARED_IMPL_HPP

#ifdef WITH_EXTRA_CHECKS
#include "function.hpp"
#endif // WITH_EXTRA_CHECKS

namespace casadi {

  template<typename Shared, typename Internal>
  void GenericShared<Shared, Internal>::count_up() {
#ifdef WITH_EXTRA_CHECKS
    casadi_assert_dev(Function::call_depth_==0);
#endif // WITH_EXTRA_CHECKS

    if (node) static_cast<Internal*>(node)->count++;

  }

  template<typename Shared, typename Internal>
  void GenericShared<Shared, Internal>::count_down() {
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

      if (--static_cast<Internal*>(node)->count == 0) {
        delete node;
        node = nullptr;
      }
    } else {
      if (--static_cast<Internal*>(node)->count == 0) {
        delete node;
        node = nullptr;
      }
    }
  }

  template<typename Shared, typename Internal>
  void GenericShared<Shared, Internal>::own(Internal* node_) {
    count_down();
    node = node_;
    count_up();
  }

  template<typename Shared, typename Internal>
  void GenericShared<Shared, Internal>::assign(Internal* node_) {
    node = node_;
  }

  template<typename Shared, typename Internal>
  std::string GenericShared<Shared, Internal>::debug_repr() const {
    if (node) {
        return node->debug_repr(node);
    } else {
        return "NULL";
    }
  }

  template<typename Shared, typename Internal>
  GenericShared<Shared, Internal>&
  GenericShared<Shared, Internal>::operator=(const GenericShared& ref) {
    // quick return if the old and new pointers point to the same object
    if (node == ref.node) return *this;

    // decrease the counter and delete if this was the last pointer
    count_down();

    // save the new pointer
    node = ref.node;
    count_up();
    return *this;
  }

  template<typename Shared, typename Internal>
  Internal* GenericShared<Shared, Internal>::get() const {
    return node;
  }

  template<typename Shared, typename Internal>
  bool GenericShared<Shared, Internal>::is_null() const {
    return node==nullptr;
  }

  template<typename Shared, typename Internal>
  Internal* GenericShared<Shared, Internal>::operator->() const {
    casadi_assert_dev(!is_null());
    return node;
  }

  template<typename Shared, typename Internal>
  void GenericShared<Shared, Internal>::swap(GenericShared& other) {
    GenericShared<Shared, Internal> temp = *this;
    *this = other;
    other = temp;
  }

  template<typename Shared, typename Internal>
  casadi_int GenericShared<Shared, Internal>::getCount() const {
    return (*this)->getCount();
  }

  template<typename Shared, typename Internal>
  GenericWeakRef<Shared, Internal>* GenericShared<Shared, Internal>::weak() {
    return (*this)->weak();
  }

  template<typename Shared, typename Internal>
  casadi_int GenericShared<Shared, Internal>::__hash__() const {
    return reinterpret_cast<casadi_int>(get());
  }

  template<typename Shared, typename Internal>
  GenericWeakRef<Shared, Internal>::GenericWeakRef(int dummy) {
    casadi_assert_dev(dummy==0);
  }

  template<typename Shared, typename Internal>
  bool GenericWeakRef<Shared, Internal>::alive() const {
    return !is_null() && (*this)->raw_ != nullptr;
  }

  template<typename Shared, typename Internal>
  Shared GenericWeakRef<Shared, Internal>::shared() const {
    Shared ret;
    if (alive()) {
      ret.own((*this)->raw_);
    }
    return ret;
  }

  template<typename Shared, typename Internal>
  bool GenericWeakRef<Shared, Internal>::shared_if_alive(Shared& shared) const {
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

  template<typename Shared, typename Internal>
  const GenericWeakRefInternal<Shared, Internal>*
  GenericWeakRef<Shared, Internal>::operator->() const {
    return static_cast<const GenericWeakRefInternal<Shared, Internal>*>(
      GenericShared<Shared, Internal>::operator->());
  }

  template<typename Shared, typename Internal>
  GenericWeakRefInternal<Shared, Internal>*
  GenericWeakRef<Shared, Internal>::operator->() {
    return static_cast<GenericWeakRefInternal<Shared, Internal>*>(
      GenericShared<Shared, Internal>::operator->());
  }

  template<typename Shared, typename Internal>
  GenericWeakRef<Shared, Internal>::GenericWeakRef(Shared shared) {
    this->own(shared.weak()->get());
  }

  template<typename Shared, typename Internal>
  GenericWeakRef<Shared, Internal>::GenericWeakRef(Internal* raw) {
    this->own(new typename Internal::weak_ref_type(raw));
  }

  template<typename Shared, typename Internal>
  void GenericWeakRef<Shared, Internal>::kill() {
    casadi_assert_dev((*this)->raw_);
    (*this)->raw_ = nullptr;
  }

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  template<typename Shared, typename Internal>
  std::shared_ptr<std::mutex> GenericWeakRef<Shared, Internal>::get_mutex() const {
    return (*this)->mutex_;
  }
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

} // namespace casadi


#endif // CASADI_GENERIC_SHARED_IMPL_HPP
