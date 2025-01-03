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


#ifndef CASADI_GENERIC_SHARED_INTERNAL_HPP
#define CASADI_GENERIC_SHARED_INTERNAL_HPP

#include "generic_shared.hpp"
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
#include <memory>
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

#ifdef CASADI_WITH_THREAD
#include <atomic>
#endif // CASADI_WITH_THREAD

namespace casadi {

  /// \cond INTERNAL
  /// Internal class for the reference counting framework, see comments on the public class.
  template<typename Shared, typename Internal>
  class GenericSharedInternal {
    friend class GenericShared<Shared, Internal>;
  public:

    /// Default constructor
    GenericSharedInternal();

    /// Copy constructor
    GenericSharedInternal(const GenericSharedInternal& node);

    /// Assignment operator
    GenericSharedInternal& operator=(const GenericSharedInternal& node);

    /// Destructor
    virtual ~GenericSharedInternal() = 0;

    /// Get the reference count
    casadi_int getCount() const;

    std::string debug_repr(const Internal*) const;

    /** \brief Get a weak reference to the object

        \identifier{1ai} */
    GenericWeakRef<Shared, Internal>* weak();

  protected:
    /** Called in the constructor of singletons to avoid that the counter reaches zero */
    void initSingleton() {
      casadi_assert_dev(static_cast<Internal*>(this)->count==0);
      static_cast<Internal*>(this)->count++;
    }

    /** Called in the destructor of singletons */
    void destroySingleton() {
      static_cast<Internal*>(this)->count--;
    }

    /// Get a shared object from the current internal object
    template<class B>
    B shared_from_this() {
      casadi_assert_dev(B::test_cast(static_cast<Internal*>(this)));
      B ret;
      ret.own(static_cast<Internal*>(this));
      return ret;
    }

    /// Get a shared object from the current internal object
    template<class B>
    const B shared_from_this() const {
      casadi_assert_dev(B::test_cast(static_cast<const Internal*>(this)));
      B ret;
      ret.own(const_cast<Internal*>(static_cast<const Internal*>(this)));
      return ret;
    }

  private:
    /// Weak pointer (non-owning) object for the object
    GenericWeakRef<Shared, Internal>* weak_ref_;
  };

  template<typename Shared, typename Internal>
  class GenericWeakRefInternal : public Internal {
  public:
    // Constructor
    GenericWeakRefInternal(Internal* raw);

    // Destructor
    ~GenericWeakRefInternal() override;

    // Raw pointer to the cached object
    Internal* raw_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    mutable std::shared_ptr<std::mutex> mutex_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
  };


  template<class A>
  A getcopy(const A& a,
             std::map<typename A::base_type*,
             typename A::internal_base_type> & already_copied) {
    A ret;
    if (!a.is_null()) {
      auto it =
          already_copied.find(const_cast<typename A::base_type*>(a.get()));
      if (it!=already_copied.end()) {
        ret.own(it->second.get());
      }
    }
    return ret;
  }

  /// \endcond


  template<typename Shared, typename Internal>
  GenericSharedInternal<Shared, Internal>::
  GenericSharedInternal(const GenericSharedInternal& node) {
    static_cast<Internal*>(this)->count = 0; // reference counter is _not_ copied
    weak_ref_ = nullptr; // nor will they have the same weak references
  }

  template<typename Shared, typename Internal>
  GenericSharedInternal<Shared, Internal>&
  GenericSharedInternal<Shared, Internal>::
  operator=(const GenericSharedInternal<Shared, Internal>& node) {
    // do _not_ copy the reference counter
    return *this;
  }

  template<typename Shared, typename Internal>
  GenericSharedInternal<Shared, Internal>::GenericSharedInternal() {
    static_cast<Internal*>(this)->count = 0;
    weak_ref_ = nullptr;
  }

  template<typename Shared, typename Internal>
  std::string GenericSharedInternal<Shared, Internal>::debug_repr(const Internal* i) const {
    // Note: i != this because of something something multiple inheritance
    return str( (casadi_int)(i)) + "/" + static_cast<const Internal*>(this)->class_name();
  }

  template<typename Shared, typename Internal>
  GenericSharedInternal<Shared, Internal>::~GenericSharedInternal() {
    #ifdef WITH_REFCOUNT_WARNINGS
    if (static_cast<Internal*>(this)->count!=0) {
      // Note that casadi_assert_warning cannot be used in destructors
      std::cerr << "Reference counting failure." <<
                   "Possible cause: Circular dependency in user code." << std::endl;
    }
    #endif // WITH_REFCOUNT_WARNINGS
    if (weak_ref_!=nullptr) {
      // Assumption: no other GenericSharedInternal instances
      // point to the same WeakRefInternal through weak_ref_
      weak_ref_->kill();
      delete weak_ref_;
      weak_ref_ = nullptr;
    }
  }

  template<typename Shared, typename Internal>
  casadi_int GenericSharedInternal<Shared, Internal>::getCount() const {
    return static_cast<const Internal*>(this)->count;
  }

  template<typename Shared, typename Internal>
  GenericWeakRef<Shared, Internal>* GenericSharedInternal<Shared, Internal>::weak() {
    if (weak_ref_==nullptr) {
      weak_ref_ = new GenericWeakRef<Shared, Internal>(static_cast<Internal*>(this));
    }
    return weak_ref_;
  }

  template<typename Shared, typename Internal>
  GenericWeakRefInternal<Shared, Internal>::GenericWeakRefInternal(Internal* raw) :
    raw_(raw)
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    , mutex_(std::make_shared<std::mutex>())
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
    {
  }

  template<typename Shared, typename Internal>
  GenericWeakRefInternal<Shared, Internal>::~GenericWeakRefInternal() {
  }


} // namespace casadi


#endif // CASADI_GENERIC_SHARED_INTERNAL_HPP
