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


#ifndef CASADI_GENERIC_SHARED_HPP
#define CASADI_GENERIC_SHARED_HPP

#include "casadi_common.hpp"
#include "exception.hpp"
#include <unordered_map>
#include <vector>
#include <cstdint>
#ifdef CASADI_WITH_THREAD
#ifdef CASADI_WITH_THREAD_MINGW
#include <mingw.mutex.h>
#else // CASADI_WITH_THREAD_MINGW
#include <mutex>
#endif // CASADI_WITH_THREAD_MINGW
#endif //CASADI_WITH_THREAD

#include <memory>

namespace casadi {

  // Forward declaration of weak reference class
  template<typename Shared, typename Internal>
  class GenericWeakRef;

  /// \cond INTERNAL
  // Forward declaration of internal classes
  template<typename Shared, typename Internal>
  class GenericSharedInternal;

  template<typename Shared, typename Internal>
  class GenericWeakRefInternal;
  /// \endcond

  template<typename Shared, typename Internal>
  class CASADI_EXPORT GenericShared {
#ifndef SWIG
    template<class B, class S> friend B shared_cast(S& A);
    template<class B, class S> friend const B shared_cast(const S& A);
#endif // SWIG

  public:
#ifndef SWIG
    /// Default constructor
    GenericShared() {
        node = nullptr;
    }

    /// Copy constructor (shallow copy)
    GenericShared(const GenericShared& ref) {
        node = ref.node;
        count_up();
    }

    /// Destructor
    ~GenericShared() {
      count_down();
    }

    /// Assignment operator
    GenericShared& operator=(const GenericShared& ref);

    /// \cond INTERNAL
    /// Assign the node to a node class pointer (or null)
    void own(Internal* node);

    /** \brief Assign the node to a node class pointer without reference counting
     *
     * improper use will cause memory leaks!

        \identifier{at} */
    void assign(Internal* node);

    /// Get a const pointer to the node
    Internal* get() const;

    /// Get the reference count
    casadi_int getCount() const;

    /// Swap content with another instance
    void swap(GenericShared& other);

    /// Access a member function or object
    Internal* operator->() const;
    /// \endcond
#endif // SWIG

    std::string debug_repr() const;


    /// Is a null pointer?
    bool is_null() const;

    /** \brief Returns a number that is unique for a given Node.

     * If the Object does not point to any node, "0" is returned.

        \identifier{av} */
    casadi_int __hash__() const;

/// \cond INTERNAL
#ifndef SWIG
    /** \brief Get a weak reference to the object

        \identifier{aw} */
    GenericWeakRef<Shared, Internal>* weak();
  protected:
    void count_up(); // increase counter of the node
    void count_down(); // decrease counter of the node
  private:
    Internal *node;
#endif // SWIG
/// \endcond
  };

  template<typename Shared, typename Internal>
  class CASADI_EXPORT GenericWeakRef : public GenericShared<Shared, Internal> {
  public:
    friend class GenericSharedInternal<Shared, Internal>;

    using GenericShared<Shared, Internal>::is_null;

    /** \brief Default constructor

        \identifier{ay} */
    GenericWeakRef(int dummy=0);

    /** \brief Construct from a shared object (also implicit type conversion)

        \identifier{az} */
    GenericWeakRef(Shared shared);

    /** \brief Get a shared (owning) reference

        \identifier{b0} */
    Shared shared() const;

    /** \brief Check if alive

        \identifier{b1} */
    bool alive() const;

    /** \brief Thread-safe alternative to alive()/shared()

        \identifier{29i} */
    bool shared_if_alive(Shared& shared) const;

    /** \brief  Access functions of the node

        \identifier{b2} */
    GenericWeakRefInternal<Shared, Internal>* operator->();

    /** \brief  Const access functions of the node

        \identifier{b3} */
    const GenericWeakRefInternal<Shared, Internal>* operator->() const;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    std::shared_ptr<std::mutex> get_mutex() const;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

#ifndef SWIG
  private:
    /** \brief Construct from a shared object (internal)

        \identifier{b4} */
    explicit GenericWeakRef(Internal* raw);

    /** \brief The shared object has been deleted

        \identifier{b5} */
    void kill();
#endif // SWIG
  };

#ifndef SWIG

  /** \brief Typecast a shared object to a base class to a shared object to a derived class,

   * cf. dynamic_cast

      \identifier{b6} */
  template<class B, class S>
  B shared_cast(S& A) {

    /// Get a pointer to the node
    typename S::internal_base_type* ptr = A.get();

    /// Create a return object
    B ret;

    /// Quick return if not allowed
    if (!B::test_cast(ptr)) return ret;

    /// Assign node of B and return
    ret.own(ptr);
    return ret;
  }

  /** \brief Typecast a shared object to a base class to a shared object to a derived class,

   * cf. dynamic_cast (const)

      \identifier{b7} */
  template<class B, class S>
  const B shared_cast(const S& A) {
    S A_copy = A;
    return shared_cast<B, S>(A_copy);
  }

#endif // SWIG

/**
 * Key is stored as a regular ref
 * Value is stored as weakref
 */
template<typename K, typename T>
class CASADI_EXPORT WeakCache {
  public:
    void tocache(const K& key, const T& f, bool needs_lock=true) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      // Add to cache
      cache_.insert(std::make_pair(key, f));
      // Remove a lost reference, if any, to prevent uncontrolled growth
      for (auto it = cache_.begin(); it!=cache_.end(); ++it) {
        if (!it->second.alive()) {
          cache_.erase(it);
          break; // just one dead reference is enough
        }
      }
    }
    /* \brief Thread-safe unique caching
    * While an incache/tocache pair in multi-threaded context is safe
    * it may lead to fresh cache entries being overwritten.
    *
    * A mutex lock_guard on the scope of an incache/tocache pair
    * may lead to deadlocks.
    *
    */
    void tocache_if_missing(const K& key, T& f) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      std::lock_guard<std::mutex> lock(mtx_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      if (!incache(key, f, false)) {
        tocache(key, f, false);
      }
    }
    bool incache(const K& key, T& f, bool needs_lock=true) const {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      auto it = cache_.find(key);
      typename T::base_type temp;
      if (it!=cache_.end() && it->second.shared_if_alive(temp)) {
        f = shared_cast<T>(temp);
        return true;
      } else {
        return false;
      }
    }
    void cache(std::vector<K>& keys, std::vector<T>& entries) const {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      std::lock_guard<std::mutex> lock(mtx_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      keys.clear();
      entries.clear();
      // Add all entries that haven't been deleted
      for (auto&& cf : cache_) {
        typename T::base_type temp;
        if (cf.second.shared_if_alive(temp)) {
          keys.push_back(cf.first);
          entries.push_back(shared_cast<T>(temp));
        }
      }
    }
  private:
    std::unordered_map<K,
      GenericWeakRef<typename T::base_type, typename T::internal_base_type>
    > cache_;
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    mutable std::mutex mtx_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
};

// gcc-12 is overzealous about use-after-free warnings
// <12 or >12 works fine
#pragma GCC diagnostic push
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ == 12)
#pragma GCC diagnostic ignored "-Wuse-after-free"
#endif

/**
 * Key is stored as a weakref
 * Value is stored as regular ref
 */
template<typename K, typename T>
class CASADI_EXPORT RevWeakCache {
  public:
    void tocache(const K& key, const T& f, bool needs_lock=true) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      // Add to cache
      const void* k = key.get();
      pre_cache_.insert(std::make_pair(k, key));
      cache_.insert(std::make_pair(k, f));
      // Remove a lost reference, if any, to prevent uncontrolled growth
      for (auto it = pre_cache_.begin(); it!=pre_cache_.end(); ++it) {
        if (!it->second.alive()) {
          pre_cache_.erase(it);
          cache_.erase(it->first);
          break; // just one dead reference is enough
        }
      }
    }
    /* \brief Thread-safe unique caching
    * While an incache/tocache pair in multi-threaded context is safe
    * it may lead to fresh cache entries being overwritten.
    *
    * A mutex lock_guard on the scope of an incache/tocache pair
    * may lead to deadlocks.
    *
    */
    void tocache_if_missing(const K& key, T& f) {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      std::lock_guard<std::mutex> lock(mtx_);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      if (!incache(key, f, false)) {
        tocache(key, f, false);
      }
    }
    bool incache(const K& key, T& f, bool needs_lock=true) const {
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
      // Safe access to cache_
      casadi::conditional_lock_guard<std::mutex> lock(mtx_, needs_lock);
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
      const void* k = key.get();
      auto it = pre_cache_.find(k);
      K temp;
      if (it!=pre_cache_.end() && it->second.shared_if_alive(temp)) {
        auto it2 = cache_.find(k);
        f = it2->second;
        return true;
      } else {
        return false;
      }
    }
  private:
    std::unordered_map<const void*,
    GenericWeakRef<typename K::base_type, typename K::internal_base_type>
    > pre_cache_;
    std::unordered_map<const void*, T> cache_;
#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    mutable std::mutex mtx_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS
};

#pragma GCC diagnostic pop

} // namespace casadi


#endif // CASADI_GENERIC_SHARED_HPP
